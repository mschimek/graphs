#include "interface.hpp"

#include <algorithm>

#include "AmsSort/AmsSort.hpp"
#include "rmat-fork/rmat/generators/select.hpp"
#include "rmat-fork/rmat/rmat.hpp"

#include "include/utils.hpp"
#include "weight_generators.hpp"

namespace graphs {
WEdgeList get_rmat_edges(const RMatParams& params,
                         WeightGeneratorConfig<Weight> wgen_config,
                         MPIComm comm) {
  using RNG = rmat::generators::select_t;
  using RMat = rmat::rmat<false>;

  const std::size_t num_threads = 1;

  std::vector<WEdgeList> edge_lists(num_threads);
  std::vector<std::size_t> edge_sizes(num_threads);

  const std::size_t num_sub_tasks = 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < num_sub_tasks; ++i) {
    const std::size_t seed = (comm.rank * (num_sub_tasks + 100)) + 1 + i;
    const std::size_t seed_scramble = seed + 1000;
    const std::size_t n = 1 << params.log_n;
    const std::size_t m = (1 << params.log_m) / num_threads;

    auto& edges = edge_lists[i];
    RNG gen(seed);
    RNG gen_scramble(seed_scramble);
    RMat r(gen_scramble, params.log_n, params.a, params.b, params.c);
    const std::size_t depth = 9u < params.log_n ? 9u : params.log_n;
    r.init(depth);
    edges.reserve(m);
    wgen_config.seed = 2 * seed;
    WeightGenerator<VId, Weight, WEdge> weight_generator{wgen_config};
    auto process_edge = [&](RMat::node src, RMat::node dst) {
      if (src == dst)
        return;
      const Weight w = weight_generator(src, dst);
      edges.emplace_back(static_cast<VId>(src), static_cast<VId>(dst), w);
      edges.emplace_back(static_cast<VId>(dst), static_cast<VId>(src), w);
    };
    r.get_edges(process_edge, m, gen);
    edge_sizes[i] = edges.size();
  }
  const std::size_t num_edges =
      std::accumulate(edge_sizes.begin(), edge_sizes.end(), 0ull);
  std::exclusive_scan(edge_sizes.begin(), edge_sizes.end(), edge_sizes.begin(),
                      0ull);
  std::vector<WEdge> edges(num_edges);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < num_sub_tasks; ++i) {
    std::size_t cur_pos = edge_sizes[i];
    auto& edge_slice = edge_lists[i];
    for (std::size_t j = 0; j < edge_slice.size(); ++j) {
      edges[cur_pos] = edge_slice[j];
      ++cur_pos;
    }
  }
  return edges;
}

template <typename T> auto get_edge_ranges(T max_edge, MPIComm comm) {
  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(T), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  std::vector<T> max_edges(comm.size);
  MPI_Allgather(&max_edge, 1, mpi_edge_type, max_edges.data(), 1, mpi_edge_type,
                comm.comm);
  MPI_Type_free(&mpi_edge_type);
  return max_edges;
}

auto remove_remaining_duplicates(WEdgeList&& edges, MPIComm comm) {
  auto max_edges = get_edge_ranges(edges.back(), comm);
  const auto org_edge_size = edges.size();
  if (comm.rank > 0) {
    const auto last_edge_prev_rank = max_edges[comm.rank - 1];
    const auto is_src_dst_equal = [&](const auto& edge) {
      return edge.get_src() == last_edge_prev_rank.get_src() &&
             edge.get_dst() == last_edge_prev_rank.get_dst();
    };
    const auto it =
        std::remove_if(edges.begin(), edges.end(), is_src_dst_equal);
    edges.erase(it, edges.end());
  }
  std::size_t num_edges = edges.size();
  if (num_edges != org_edge_size) {
    std::cout << (org_edge_size - num_edges) << std::endl;
  }
  kagen::KaGenResult res;
  res.edges.resize(num_edges);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    res.edges[i] = {edges[i].get_src(), edges[i].get_dst()};
  }
  res.vertex_range.first = edges.front().get_src();
  res.vertex_range.second = edges.back().get_src();
  // final verification that there are no duplicates:
  auto SrcDstSort = [](const auto& lhs, const auto& rhs) {
    return lhs < rhs;
  };
  auto SrcDstEqual = [](const auto& lhs, const auto& rhs) {
    return lhs == rhs;
  };
  // verification
  if (!std::is_sorted(res.edges.begin(), res.edges.end(), SrcDstSort)) {
    std::cout << "edge are not sorted after remaining duplicate removal!"
              << std::endl;
    std::abort();
  }
  if (std::adjacent_find(res.edges.begin(), res.edges.end()) != res.edges.end()) {
    std::cout << "edge are not unique after remaining duplicate removal!"
              << std::endl;
    std::abort();
  }
  const auto new_max_edges = get_edge_ranges(res.edges.back(), comm);
  if (comm.rank > 0) {
    const auto max_edge_prev = new_max_edges[comm.rank - 1];
    if (std::get<0>(max_edge_prev) == std::get<0>(res.edges.front()) &&
        std::get<1>(max_edge_prev) == std::get<1>(res.edges.front())) {
      std::cout
          << "edge are not globally sorted after remaining duplicate removal!"
          << std::endl;
      std::abort();
    }
  }

  return res;
}

std::pair<WEdgeList, VertexRange>
get_rmat_edges_evenly_distributed(const RMatParams& params,
                                  WeightGeneratorConfig<Weight> wgen_config,
                                  bool remove_duplicates, MPIComm comm) {
  auto edges = get_rmat_edges(params, wgen_config, comm);
  // edges = mpi::par_sort(std::move(edges), ctx, LexicOrder{});
  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(WEdge), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  std::mt19937_64 gen(comm.rank * 100);
  auto SrcDstSort = [](const WEdge& lhs, const WEdge& rhs) {
    return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) <
           std::make_tuple(rhs.get_src(), rhs.get_dst(), lhs.get_weight());
  };
  const int num_levels = 2;
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstSort);
  if (remove_duplicates) {
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  }

  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstSort);
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  MPI_Type_free(&mpi_edge_type);


  const VId v_min = !edges.empty() ? edges.front().get_src() : -1;
  const VId v_max = !edges.empty() ? edges.back().get_src() : -1;
  auto res = remove_remaining_duplicates(std::move(edges), comm);
  UniformRandomWeightGenerator<VId, Weight, WEdge> w_gen{wgen_config};
  return add_weights_rmat(w_gen, std::move(res), comm);
}
} // namespace graphs
