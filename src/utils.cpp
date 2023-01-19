#include "include/utils.hpp"
#include "interface.hpp"
#include "ips4o/ips4o.hpp"
#include "memory_utils.hpp"
#include <numeric>
#include <sstream>

namespace graphs {
void remove_upside_down(std::vector<WEdge>& edges, const VertexRange& range) {
  auto it = std::remove_if(edges.begin(), edges.end(), [&](const WEdge& edge) {
    return edge.get_src() < range.first || edge.get_src() > range.second;
  });
  edges.erase(it, edges.end());
  edges.shrink_to_fit();
}

namespace internal {
class MessageBuffers {
public:
  MessageBuffers(MPIComm comm) : comm_{comm}, data_for_pe_(comm.size) {}
  void add(const WEdge& elem, int pe) { data_for_pe_[pe].push_back(elem); }
  std::vector<WEdge> exchange() {
    const int size = data_for_pe_.size();

    std::vector<WEdge> send_buf;
    std::vector<WEdge> recv_buf;
    std::vector<int> send_counts(size);
    std::vector<int> recv_counts(size);
    std::vector<int> send_displs(size);
    std::vector<int> recv_displs(size);
    for (size_t i = 0; i < send_counts.size(); ++i) {
      send_counts[i] = data_for_pe_[i].size();
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(),
                        send_displs.begin(), 0);
    const std::size_t total_send_count =
        send_displs.back() + send_counts.back();
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                 comm_.comm);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                        recv_displs.begin(), 0);
    const std::size_t total_recv_count =
        recv_displs.back() + recv_counts.back();

    send_buf.reserve(total_send_count);
    memory_stats().print("in exchange");
    for (std::size_t i = 0; i < send_counts.size(); ++i) {
      for (const auto& elem : data_for_pe_[i]) {
        send_buf.push_back(elem);
      }
      data_for_pe_[i].clear();
      data_for_pe_[i].shrink_to_fit();
    }

    recv_buf.resize(total_recv_count);
    const auto wedge_mpi_type = WEdge::MPI_Type{};
    MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(),
                  wedge_mpi_type.get_mpi_type(), recv_buf.data(),
                  recv_counts.data(), recv_displs.data(),
                  wedge_mpi_type.get_mpi_type(), comm_.comm);
    return recv_buf;
  }

private:
  MPIComm comm_;
  std::vector<std::vector<WEdge>> data_for_pe_;
};

std::vector<VertexRange> get_ranges(const VertexRange local_range,
                                    MPIComm comm) {
  std::vector<VertexRange> ranges(comm.size);
  ranges[comm.rank] = local_range;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(),
                sizeof(VertexRange), MPI_BYTE, comm.comm);
  // if (comm.rank == 0) {
  //   for (const auto [first, second] : ranges) {
  //     std::cout << first << ", " << second << std::endl;
  //   }
  // }
  return ranges;
}

std::vector<EdgeRange> get_edge_ranges(const EdgeRange local_range,
                                       MPIComm comm) {
  std::vector<EdgeRange> ranges(comm.size);
  ranges[comm.rank] = local_range;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(),
                sizeof(EdgeRange), MPI_BYTE, comm.comm);
  return ranges;
}

int get_home_pe(VId v, const std::vector<VertexRange>& ranges) {
  const auto is_v_smaller_than_range =
      [](const VId v, const VertexRange& range) { return v < range.first; };
  auto it = std::upper_bound(ranges.begin(), ranges.end(), v,
                             is_v_smaller_than_range);
  if (it == ranges.begin()) {
    std::cerr << "invalid home pe query" << std::endl;
    std::abort();
  }
  return std::distance(ranges.begin(), it) - 1;
}

int get_home_pe(Edge e, const std::vector<EdgeRange>& ranges) {
  const auto is_e_smaller_than_range = [](const Edge& e,
                                          const EdgeRange& range) {
    return std::make_pair(e.get_src(), e.get_dst()) <
           std::make_pair(range.first.get_src(), range.first.get_dst());
  };
  auto it = std::upper_bound(ranges.begin(), ranges.end(), e,
                             is_e_smaller_than_range);
  if (it == ranges.begin()) {
    std::cerr << "invalid home pe query" << std::endl;
    std::abort();
  }
  return std::distance(ranges.begin(), it) - 1;
}

WEdgeList
get_remote_edges_pointing_to_pe(const WEdgeList& edges,
                                const std::vector<VertexRange>& ranges,
                                MPIComm comm) {
  MessageBuffers msg_buffers(comm);
  for (const auto& edge : edges) {
    int pe = get_home_pe(edge.get_dst(), ranges);
    if (pe == comm.rank)
      continue;
    msg_buffers.add(edge, pe);
  }
  return msg_buffers.exchange();
}

WEdgeList get_remote_edges_pointing_to_pe_pseudo_inplace(
    WEdgeList& edges, const std::vector<VertexRange>& ranges, MPIComm comm) {

  ips4o::parallel::sort(edges.begin(), edges.end(), DstSrcOrder<WEdge>{});
  std::vector<int> send_counts(comm.size);
  for (const auto& edge : edges) {
    int pe = get_home_pe(edge.get_dst(), ranges);
    ++send_counts[pe];
  }
  std::vector<int> recv_counts(comm.size);
  std::vector<int> send_displs(comm.size);
  std::vector<int> recv_displs(comm.size);

  std::exclusive_scan(send_counts.begin(), send_counts.end(),
                      send_displs.begin(), 0);
  const std::size_t total_send_count = send_displs.back() + send_counts.back();
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               comm.comm);
  std::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                      recv_displs.begin(), 0);
  const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();

  WEdgeList recv_buf(total_recv_count);
  memory_stats().print("in exchange");
  const auto wedge_mpi_type = WEdge::MPI_Type{};
  MPI_Alltoallv(edges.data(), send_counts.data(), send_displs.data(),
                wedge_mpi_type.get_mpi_type(), recv_buf.data(),
                recv_counts.data(), recv_displs.data(),
                wedge_mpi_type.get_mpi_type(), comm.comm);
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstOrder<WEdge>{});
  return recv_buf;
}

WEdgeList get_remote_edges_pointing_to_pe_pseudo_inplace(
    WEdgeList& edges, const std::vector<EdgeRange>& ranges, MPIComm comm) {

  ips4o::parallel::sort(edges.begin(), edges.end(), DstSrcOrder<WEdge>{});
  std::vector<int> send_counts(comm.size);
  for (const auto& edge : edges) {
    int pe = get_home_pe(Edge{edge.get_dst(), edge.get_src()}, ranges);
    ++send_counts[pe];
  }
  std::vector<int> recv_counts(comm.size);
  std::vector<int> send_displs(comm.size);
  std::vector<int> recv_displs(comm.size);

  std::exclusive_scan(send_counts.begin(), send_counts.end(),
                      send_displs.begin(), 0);
  const std::size_t total_send_count = send_displs.back() + send_counts.back();
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               comm.comm);
  std::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                      recv_displs.begin(), 0);
  const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();

  WEdgeList recv_buf(total_recv_count);
  memory_stats().print("in exchange");
  const auto wedge_mpi_type = WEdge::MPI_Type{};
  MPI_Alltoallv(edges.data(), send_counts.data(), send_displs.data(),
                wedge_mpi_type.get_mpi_type(), recv_buf.data(),
                recv_counts.data(), recv_displs.data(),
                wedge_mpi_type.get_mpi_type(), comm.comm);
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstOrder<WEdge>{});
  return recv_buf;
}

std::size_t remove_duplicates(WEdgeList& edges) {
  auto it = std::unique(edges.begin(), edges.end(), SrcDstWeightEqual<WEdge>{});
  std::size_t num_duplicates = std::distance(it, edges.end());
  edges.erase(it, edges.end());
  return num_duplicates;
}

std::size_t add_missing_edges(WEdgeList& edges, const WEdgeList& remote_edges) {
  std::vector<WEdge> missing_edges;
  const auto comp = SrcDstOrder<WEdge>{};
  for (const auto& remote_edge : remote_edges) {
    const WEdge flipped_edge{remote_edge.get_dst(), remote_edge.get_src(),
                             remote_edge.get_weight()};
    const auto it =
        std::lower_bound(edges.begin(), edges.end(), flipped_edge, comp);
    if (it == edges.end()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->get_src() != flipped_edge.get_src() ||
               it->get_dst() != flipped_edge.get_dst()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->weight != remote_edge.weight) {
      it->weight = comp(*it, remote_edge) ? it->weight : remote_edge.weight;
    }
  }
  if (!missing_edges.empty()) {
    edges.insert(edges.end(), missing_edges.begin(), missing_edges.end());
    std::sort(edges.begin(), edges.end(), SrcDstWeightOrder<WEdge>{});
  }
  return missing_edges.size();
}

// assertion no duplicates in edges
std::size_t add_missing_edges_via_merging(WEdgeList& edges,
                                          const WEdgeList& remote_edges) {
  std::vector<WEdge> missing_edges;
  const auto comp = SrcDstOrder<WEdge>{};
  auto it = edges.begin();
  MPIComm comm;
  for (std::size_t i = 0; i < remote_edges.size(); ++i) {
    const auto& remote_edge = remote_edges[i];
    const WEdge flipped_edge{remote_edge.get_dst(), remote_edge.get_src(),
                             remote_edge.get_weight()};
    for (; it != edges.end() && comp(*it, flipped_edge); ++it)
      ;
    // now it points to the first element of edges which is equal to or greater
    // than flipped_edge
    if (it == edges.end()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->get_src() != flipped_edge.get_src() ||
               it->get_dst() != flipped_edge.get_dst()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->weight != remote_edge.weight) {
      it->weight = comp(*it, remote_edge) ? it->weight : remote_edge.weight;
    }
  }
  if (!missing_edges.empty()) {
    edges.insert(edges.end(), missing_edges.begin(), missing_edges.end());
    ips4o::parallel::sort(edges.begin(), edges.end(),
                          SrcDstWeightOrder<WEdge>{});
  }
  return missing_edges.size();
}
bool is_local(VId v, const VertexRange& range) {
  return range.first <= v && v <= range.second;
}

bool is_local(WEdge edge, const VertexRange& range) {
  return is_local(edge.get_src(), range) && is_local(edge.get_dst(), range);
}

std::size_t add_missing_local_edges(WEdgeList& edges,
                                    const VertexRange& local_range) {
  std::vector<WEdge> missing_edges;
  const auto comp = SrcDstOrder<WEdge>{};
  for (auto& edge : edges) {
    if (!is_local(edge, local_range)) {
      continue;
    }
    const WEdge flipped_edge{edge.get_dst(), edge.get_src(), edge.get_weight()};
    const auto it =
        std::lower_bound(edges.begin(), edges.end(), flipped_edge, comp);
    if (it == edges.end()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->get_src() != flipped_edge.get_src() ||
               it->get_dst() != flipped_edge.get_dst()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->weight != edge.weight) {
      edge.weight = comp(edge, *it) ? edge.weight : it->weight;
    }
  }
  if (!missing_edges.empty()) {
    edges.insert(edges.end(), missing_edges.begin(), missing_edges.end());
    std::sort(edges.begin(), edges.end(), SrcDstWeightOrder<WEdge>{});
  }
  return missing_edges.size();
}

std::uint64_t allreduce_sum(const std::uint64_t& send_elem, MPIComm comm) {
  std::uint64_t recv_elem;
  MPI_Allreduce(&send_elem, &recv_elem, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  return recv_elem;
}

} // namespace internal

// assumption: Vertices do not span multiple PEs
// assumption: Edges are already (locally) sorted
void repair_edges(WEdgeList& edges, const VertexRange& local_range,
                  MPIComm comm) {
  using namespace internal;
  memory_stats().print("before repair");
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstOrder<WEdge>{});

  const auto ranges = get_ranges(local_range, comm);

  // if(comm.rank == 0) {
  //   for(const auto [vmin, vmax] : ranges) {
  //     std::cout << vmin << ", " << vmax << std::endl;
  //   }
  // }
  const std::uint64_t num_duplicates = remove_duplicates(edges);
  auto remote_edges =
      get_remote_edges_pointing_to_pe_pseudo_inplace(edges, ranges, comm);
  memory_stats().print("after exchanging edges");
  ips4o::parallel::sort(remote_edges.begin(), remote_edges.end(),
                        DstSrcOrder<WEdge>{});
  memory_stats().print("after exchanging edges I");
  memory_stats().print("after exchanging edges II");
  const std::uint64_t num_missing_local_edges =
      add_missing_local_edges(edges, local_range);
  memory_stats().print("after exchanging edges III");
  // auto copy = edges;
  const std::uint64_t num_missing_edges =
      add_missing_edges_via_merging(edges, remote_edges);
  // const std::uint64_t num_missing_edges_ =
  //     add_missing_edges(copy, remote_edges);
  // if(copy != edges || num_missing_edges != num_missing_edges_) {
  //   std::abort();
  // }
  remove_duplicates(edges);
  memory_stats().print("after exchanging edges IV");

  std::stringstream sstream;
  sstream << "rank: " << comm.rank
          << " num_duplicates: " << allreduce_sum(num_duplicates, comm)
          << " num_missing_local_edges: "
          << allreduce_sum(num_missing_local_edges, comm)
          << " num_missing_edges: " << allreduce_sum(num_missing_edges, comm)
          << " remote_edges: " << allreduce_sum(remote_edges.size(), comm)
          << std::endl;
  if (comm.rank == 0) {
    std::cout << sstream.str() << std::endl;
  }
}

std::size_t get_next_pow_two_with_exp_divisible_by(std::size_t i,
                                                   std::size_t divisor) {
  if (i <= 1) {
    return 1ull << divisor;
  }
  const std::size_t exponent = static_cast<std::size_t>(std::log2(i));
  const std::size_t remainder = exponent % divisor;
  const std::size_t offset = remainder == 0 ? 0 : (divisor - remainder);
  const std::size_t exponent_divisible_by = exponent + offset;
  return 1ull << exponent_divisible_by;
}

std::vector<WEdge> get_local_weighted_edges(
    UniformRandomWeightGenerator<VId, Weight, WEdge>& w_gen,
    const kagen::KaGenResult& kagen_result) {
  const auto& [edges, vertex_range] = kagen_result;
  std::vector<WEdge> w_edges;
  w_edges.reserve(edges.size());
  for (const auto& [src, dst] : edges) {
    w_edges.emplace_back(src, dst, w_gen(src, dst));
  }
  return w_edges;
}
void set_global_weighted_edges(WEdgeList& remote_edges, WEdgeList& own_edges, MPIComm comm) {
  ips4o::parallel::sort(remote_edges.begin(), remote_edges.end(),
                        DstSrcOrder<WEdge>{});

  //MPI_Comm comm;
  //execute_in_order(comm, [&](){
  if (remote_edges.size() != own_edges.size()) {
    std::cout << remote_edges.size() << " " << own_edges.size() << std::endl;
    std::cout << "error: -----" << std::endl;
    for(const auto& edge_ : remote_edges) {
      std::cout << edge_ << std::endl;
      const auto it = std::find_if(own_edges.begin(), own_edges.end(), [&](const auto& edge) {
              return (edge.get_src() == edge_.get_dst() && edge.get_dst() == edge_.get_src());
          });
      if(it == own_edges.end()) {
        std::cout << "here" << std::endl;
      }
    }
    std::cout << "-----" << std::endl;
    for(const auto& edges : own_edges) {
      std::cout << edges << std::endl;
    }
    std::cout << "edges do not have equal size!" << std::endl;
    std::abort();
  }
  //});
  auto comp = SrcDstOrder<WEdge>{};
  for (std::size_t i = 0; i < own_edges.size(); ++i) {
    auto& edge = own_edges[i];
    const auto& remote_edge = remote_edges[i];
    if (edge.get_src() != remote_edge.get_dst() ||
        edge.get_dst() != remote_edge.get_src()) {
      std::cout << "edges are corrupt!" << std::endl;
      std::abort();
    }
    const auto weight =
        comp(edge, remote_edge) ? edge.get_weight() : remote_edge.get_weight();
    edge.set_weight(weight);
  }
}

std::pair<std::vector<WEdge>, VertexRange>
add_weights(UniformRandomWeightGenerator<VId, Weight, WEdge>& w_gen,
            kagen::KaGenResult&& kagen_result, MPIComm comm) {
  // assert: all back edges are existant and there are no duplicates
  std::pair<std::vector<WEdge>, VertexRange> result;
  result.first = get_local_weighted_edges(w_gen, kagen_result);
  result.second = kagen_result.vertex_range;
  auto& [w_edges, local_range] = result;
  const auto ranges = internal::get_ranges(local_range, comm);

  kagen::KaGenResult{{}, {}} = std::move(kagen_result); // clear storage
  auto remote_edges = internal::get_remote_edges_pointing_to_pe_pseudo_inplace(
      w_edges, ranges, comm);
  set_global_weighted_edges(remote_edges, w_edges, comm);
  return result;
}

std::pair<std::vector<WEdge>, VertexRange>
add_weights_rmat(UniformRandomWeightGenerator<VId, Weight, WEdge>& w_gen,
                 kagen::KaGenResult&& kagen_result, MPIComm comm) {
  // assert: all back edges are existant and there are no duplicates
  std::pair<std::vector<WEdge>, VertexRange> result;
  result.first = get_local_weighted_edges(w_gen, kagen_result);
  result.second = kagen_result.vertex_range;
  auto& [w_edges, local_range] = result;
  Edge min_edge{w_edges.front().get_src(), w_edges.front().get_dst()};
  Edge max_edge{w_edges.back().get_src(), w_edges.back().get_dst()};
  EdgeRange edge_range{min_edge, max_edge};
  const auto ranges = internal::get_edge_ranges(edge_range, comm);
  //execute_in_order(comm, [&]() {
  //  if (comm.rank == 0) {
  //    for (const auto range : ranges) {
  //      std::cout << range.first << " " << range.second << std::endl;
  //    }
  //  }
  //  for (const auto& edge : kagen_result.edges) {
  //    std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << std::endl;
  //  }
  //});
  kagen::KaGenResult{{}, {}} = std::move(kagen_result); // clear storage
  auto remote_edges = internal::get_remote_edges_pointing_to_pe_pseudo_inplace(
      w_edges, ranges, comm);
  set_global_weighted_edges(remote_edges, w_edges, comm);
  return result;
}
} // namespace graphs
