#include "include/utils.hpp"
#include "interface.hpp"
#include <numeric>
#include <sstream>

namespace graphs {
void remove_upside_down(WEdgeList& edges, const VertexRange& range) {
  auto it = std::remove_if(edges.begin(), edges.end(), [&](const WEdge& edge) {
    return edge.src < range.first || edge.src > range.second;
  });
  edges.erase(it, edges.end());
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

WEdgeList
get_remote_edges_pointing_to_pe(const WEdgeList& edges,
                                const std::vector<VertexRange>& ranges,
                                MPIComm comm) {
  MessageBuffers msg_buffers(comm);
  for (const auto& edge : edges) {
    int pe = get_home_pe(edge.dst, ranges);
    if (pe == comm.rank)
      continue;
    msg_buffers.add(edge, pe);
  }
  return msg_buffers.exchange();
}

std::size_t remove_duplicates(WEdgeList& edges) {
  auto it = std::unique(edges.begin(), edges.end());
  std::size_t num_duplicates = std::distance(it, edges.end());
  edges.erase(it, edges.end());
  return num_duplicates;
}

std::size_t add_missing_edges(WEdgeList& edges, const WEdgeList& remote_edges) {
  std::vector<WEdge> missing_edges;
  const auto comp = SrcDstOrder{};
  for (const auto& remote_edge : remote_edges) {
    const WEdge flipped_edge{remote_edge.dst, remote_edge.src,
                             remote_edge.weight};
    const auto it = std::lower_bound(edges.begin(), edges.end(), flipped_edge,
                                     comp);
    if (it == edges.end()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->src != flipped_edge.src || it->dst != flipped_edge.dst) {
      missing_edges.push_back(flipped_edge);
    } else if (it->weight != remote_edge.weight) {
      it->weight = comp(*it, remote_edge) ? it->weight : remote_edge.weight;
    }
  }
  if (!missing_edges.empty()) {
    edges.insert(edges.end(), missing_edges.begin(), missing_edges.end());
    std::sort(edges.begin(), edges.end(), SrcDstWeightOrder{});
  }
  return missing_edges.size();
}
bool is_local(VId v, const VertexRange& range) {
  return range.first <= v && v <= range.second;
}

bool is_local(WEdge edge, const VertexRange& range) {
  return is_local(edge.src, range) && is_local(edge.dst, range);
}

std::size_t add_missing_local_edges(WEdgeList& edges,
                                    const VertexRange& local_range) {
  std::vector<WEdge> missing_edges;
  const auto comp = SrcDstOrder{};
  for (auto& edge : edges) {
    if (!is_local(edge, local_range)) {
      continue;
    }
    const WEdge flipped_edge{edge.dst, edge.src, edge.weight};
    const auto it = std::lower_bound(edges.begin(), edges.end(), flipped_edge,
                                     comp);
    if (it == edges.end()) {
      missing_edges.push_back(flipped_edge);
    } else if (it->src != flipped_edge.src || it->dst != flipped_edge.dst) {
      missing_edges.push_back(flipped_edge);
    } else if (it->weight != edge.weight) {
      edge.weight = comp(edge, *it) ? edge.weight : it->weight;
    }
  }
  if (!missing_edges.empty()) {
    edges.insert(edges.end(), missing_edges.begin(), missing_edges.end());
    std::sort(edges.begin(), edges.end(), SrcDstWeightOrder{});
  }
  return missing_edges.size();
}

std::uint64_t allreduce_sum(const std::uint64_t& send_elem, MPIComm comm) {
  std::uint64_t recv_elem;
  MPI_Allreduce(&send_elem, &recv_elem, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  return recv_elem;
}

} // namespace internal

void repair_edges(WEdgeList& edges, const VertexRange& local_range,
                  MPIComm comm) {
  using namespace internal;
  std::sort(edges.begin(), edges.end(), SrcDstOrder{});

  const auto ranges = get_ranges(local_range, comm);

  // if(comm.rank == 0) {
  //   for(const auto [vmin, vmax] : ranges) {
  //     std::cout << vmin << ", " << vmax << std::endl;
  //   }
  // }
  auto remote_edges = get_remote_edges_pointing_to_pe(edges, ranges, comm);
  std::sort(remote_edges.begin(), remote_edges.end(), SrcDstOrder{});
  const std::uint64_t num_duplicates = remove_duplicates(edges);
  const std::uint64_t num_missing_local_edges =
      add_missing_local_edges(edges, local_range);
  const std::uint64_t num_missing_edges =
      add_missing_edges(edges, remote_edges);

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
} // namespace graphs
