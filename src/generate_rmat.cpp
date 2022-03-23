#include "interface.hpp"

#include <algorithm>

#include "AmsSort/AmsSort.hpp"
#include "rmat-fork/rmat/generators/select.hpp"
#include "rmat-fork/rmat/rmat.hpp"

#include "weight_generators.hpp"

namespace graphs {
WEdgeList get_rmat_edges(const RMatParams& params, MPIComm comm) {
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
    WeightGenerator<VId, Weight, WEdge> weight_generator{};
    auto process_edge = [&](RMat::node src, RMat::node dst) {
      if (src == dst) return;
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

//// WEdgeList get_rmat_edges(const RMatParams& params, VId v_begin,
////                                   VId v_end, bool remove_duplicates) {
////   mpi::MPIContext ctx;
////   auto edges = get_rmat_edges(params);
////
////   VertexRange own_range{v_begin, v_end, ctx.rank()};
////   const std::vector<VertexRange> ranges = mpi::allgather(own_range);
////   SEQ_EX(ctx, {
////     PRINT_VECTOR(edges);
////     PRINT_VECTOR(ranges);
////   });
////   PEIDToVector<WEdge> pe_edges_send;
////   PEIDToVector<WEdge> pe_edges_recv;
////   for (const auto& edge : edges) {
////     const auto [from, to, w] = edge;
////     pe_edges_send[get_pe(from, ranges)].push_back(edge);
////     pe_edges_send[get_pe(to, ranges)].emplace_back(to, from,
////                                                    w);  // reverse edge
////   }
////   edges.clear();
////   mpi::sparse_all_to_all(pe_edges_send, pe_edges_recv, ctx);
////   pe_edges_send.clear();
////   edges = mpi::make_contiguous(pe_edges_recv);
////
////   if (remove_duplicates) {
////     std::sort(edges.begin(), edges.end(), LexicOrder{});
////     edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
////   }
////   return edges;
//// }
////
std::pair<WEdgeList, VertexRange> get_rmat_edges_evenly_distributed(
    const RMatParams& params, bool remove_duplicates, MPIComm comm) {
  auto edges = get_rmat_edges(params, comm);
  // edges = mpi::par_sort(std::move(edges), ctx, LexicOrder{});
  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(WEdge), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  std::mt19937_64 gen(comm.rank * 100);
  auto SrcDstSort = [](const WEdge& lhs, const WEdge& rhs) {
    return std::tie(lhs.src, lhs.dst, lhs.weight) <
           std::tie(rhs.src, rhs.dst, lhs.weight);
  };
  const int num_levels = 2;
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstSort);
  if (remove_duplicates) {
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  }

  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstSort);
  MPI_Type_free(&mpi_edge_type);
  const VId v_min = !edges.empty() ? edges.front().src : -1;
  const VId v_max = !edges.empty() ? edges.back().src : -1;
  return std::make_pair(std::move(edges), VertexRange{v_min, v_max});
}
}  // namespace graphs
