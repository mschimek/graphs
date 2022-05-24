#pragma once

#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>

#include "AmsSort/AmsSort.hpp"

#include "interface.hpp"

namespace graphs {

inline std::size_t compute_number_edges(const std::string& infile,
                                        std::size_t num_VId_bytes,
                                        MPIComm comm = MPIComm{}) {
  std::ifstream in{infile.c_str(), std::ios::binary | std::ios::in};
  std::ifstream end_of_file(infile.c_str(), std::ios::binary | std::ios::ate);
  std::fstream::pos_type file_size;
  if (end_of_file) {
    file_size = end_of_file.tellg();
  }
  std::size_t step_size = 2 * num_VId_bytes + 1; // src, dst, weight
  const std::size_t num_edges_global = file_size / step_size;
  const std::size_t num_edges_local = num_edges_global / comm.size;
  const std::size_t own_num_edges =
      comm.rank + 1 < comm.size
          ? num_edges_local
          : num_edges_local + num_edges_global % comm.size;
  return own_num_edges;
}

template <typename OutIt>
std::size_t read_weighted_binary(const std::string& infile, OutIt out_it,
                                 std::size_t num_VId_bytes,
                                 MPIComm comm = MPIComm{}) {
  if (comm.rank == 0) {
    std::cout << infile << std::endl;
  }
  using EdgeType = typename std::iterator_traits<OutIt>::value_type;
  std::ifstream in{infile.c_str(), std::ios::binary | std::ios::in};
  std::ifstream end_of_file(infile.c_str(), std::ios::binary | std::ios::ate);
  std::fstream::pos_type file_size;
  if (end_of_file) {
    file_size = end_of_file.tellg();
  }
  std::size_t step_size = 2 * num_VId_bytes + 1; // src, dst, weight
  const std::size_t num_edges = file_size / step_size;
  const std::size_t chunk_size = num_edges / comm.size;
  const std::size_t local_chunk_size = comm.rank + 1 < comm.size
                                           ? chunk_size
                                           : chunk_size + num_edges % comm.size;
  std::vector<char> data(local_chunk_size * step_size);
  in.seekg(0, std::ios::beg);
  in.seekg(chunk_size * step_size * comm.rank);
  in.read(reinterpret_cast<char*>(data.data()), local_chunk_size * step_size);
  // std::cout << "rank: " << comm.rank << " file size: " << file_size << " "
  //            << sizeof(std::fstream::pos_type)
  //           << " local_chunk_size: " << local_chunk_size << std::endl;
  for (std::size_t i = 0; i < data.size(); i += step_size) {
    auto it = data.data() + i;
    VId src = 0;
    VId dst = 0;
    uint8_t w = 0;
    std::memcpy(&src, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&dst, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&w, it, 1);
    EdgeType e;
    e.set_src(src);
    e.set_dst(dst);
    e.set_weight(w);
    // std::cout << e << std::endl;
    *(out_it) = e;
    ++out_it;
  }
  return local_chunk_size;
}

namespace helpers {
using EdgeRange = std::pair<Edge, Edge>;
inline std::vector<EdgeRange> get_ranges(const EdgeRange local_range,
                                         MPIComm comm) {
  std::vector<EdgeRange> ranges(comm.size);
  ranges[comm.rank] = local_range;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(),
                sizeof(VertexRange), MPI_BYTE, comm.comm);
  return ranges;
}

template <typename Container>
inline void delete_duplicate_edge(Container& edges, MPIComm comm) {
  Edge min_edge{edges.front().get_src(), edges.front().get_dst()};
  Edge max_edge{edges.back().get_src(), edges.back().get_dst()};
  EdgeRange range{min_edge, max_edge};
  const auto ranges = get_ranges(range, comm);
  if (comm.rank == 0)
    return;
  const auto prev_max_edge = ranges[comm.rank - 1].second;
  if (prev_max_edge.get_src() == edges.front().get_src() &&
      prev_max_edge.get_dst() == edges.front().get_dst()) {
    edges.erase(edges.begin());
  }
}
inline std::vector<VertexRange> get_ranges(const VertexRange local_range,
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

inline int get_home_pe(VId v, const std::vector<VertexRange>& ranges) {
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

inline WEdgeList get_remote_edges_pointing_to_pe_pseudo_inplace(
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
  const auto wedge_mpi_type = WEdge::MPI_Type{};
  MPI_Alltoallv(edges.data(), send_counts.data(), send_displs.data(),
                wedge_mpi_type.get_mpi_type(), recv_buf.data(),
                recv_counts.data(), recv_displs.data(),
                wedge_mpi_type.get_mpi_type(), comm.comm);
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstOrder<WEdge>{});
  return recv_buf;
}

inline std::size_t remove_duplicates(WEdgeList& edges) {
  auto it = std::unique(edges.begin(), edges.end(), SrcDstWeightEqual<WEdge>{});
  std::size_t num_duplicates = std::distance(it, edges.end());
  edges.erase(it, edges.end());
  return num_duplicates;
}

inline std::size_t add_missing_edges(WEdgeList& edges, const WEdgeList& remote_edges) {
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
inline std::size_t add_missing_edges_via_merging(WEdgeList& edges,
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
    ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstWeightOrder<WEdge>{});
  }
  return missing_edges.size();
}
inline bool is_local(VId v, const VertexRange& range) {
  return range.first <= v && v <= range.second;
}

inline bool is_local(WEdge edge, const VertexRange& range) {
  return is_local(edge.get_src(), range) && is_local(edge.get_dst(), range);
}

inline std::size_t add_missing_local_edges(WEdgeList& edges,
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

inline std::uint64_t allreduce_sum(const std::uint64_t& send_elem, MPIComm comm) {
  std::uint64_t recv_elem;
  MPI_Allreduce(&send_elem, &recv_elem, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  return recv_elem;
}


// assumption: Vertices do not span multiple PEs
// assumption: Edges are already (locally) sorted
inline void repair_edges(WEdgeList& edges, const VertexRange& local_range,
                  MPIComm comm) {
  using namespace helpers;
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstOrder<WEdge>{});

  const auto ranges = helpers::get_ranges(local_range, comm);

  // if(comm.rank == 0) {
  //   for(const auto [vmin, vmax] : ranges) {
  //     std::cout << vmin << ", " << vmax << std::endl;
  //   }
  // }
  const std::uint64_t num_duplicates = remove_duplicates(edges);
  auto remote_edges =
      get_remote_edges_pointing_to_pe_pseudo_inplace(edges, ranges, comm);
  ips4o::parallel::sort(remote_edges.begin(), remote_edges.end(),
                        DstSrcOrder<WEdge>{});
  const std::uint64_t num_missing_local_edges =
      add_missing_local_edges(edges, local_range);
  // auto copy = edges;
  const std::uint64_t num_missing_edges =
      add_missing_edges_via_merging(edges, remote_edges);
  // const std::uint64_t num_missing_edges_ =
  //     add_missing_edges(copy, remote_edges);
  // if(copy != edges || num_missing_edges != num_missing_edges_) {
  //   std::abort();
  // }

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
} // namespace helpers

template <typename Container, typename Infiles>
std::pair<Container, VertexRange>
read_weighted_binaries(const Infiles& infiles, std::size_t num_VId_bytes,
                       MPIComm comm = MPIComm{}) {
  using EdgeType = typename Container::value_type;
  std::size_t num_own_edges = 0;
  for (const auto& infile : infiles) {
    num_own_edges += compute_number_edges(infile, num_VId_bytes, comm);
  }
  std::pair<Container, VertexRange> res;
  Container& edges = res.first;
  VertexRange& vertex_range = res.second;
  edges.resize(num_own_edges);
  auto it = edges.begin();
  for (const auto& infile : infiles) {
    const std::size_t num_read_edges =
        read_weighted_binary(infile, it, num_VId_bytes, comm);
    it += num_read_edges;
  }
  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(EdgeType), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  std::mt19937_64 gen(comm.rank * 100);
  auto SrcDstWeightOrder = [](const EdgeType& lhs, const EdgeType& rhs) {
    return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) <
           std::make_tuple(rhs.get_src(), rhs.get_dst(), rhs.get_weight());
  };
  auto SrcDstEqual = [](const EdgeType& lhs, const EdgeType& rhs) {
    return std::make_tuple(lhs.get_src(), lhs.get_dst()) ==
           std::make_tuple(rhs.get_src(), rhs.get_dst());
  };
  const int num_levels = 2;
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm,
                 SrcDstWeightOrder);
  edges.erase(std::unique(edges.begin(), edges.end(), SrcDstEqual),
              edges.end());
  helpers::delete_duplicate_edge(edges, comm);
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm,
                 SrcDstWeightOrder);
  if (edges.empty()) {
    vertex_range.first = -1;
    vertex_range.second = -1;
  } else {
    vertex_range.first = edges.front().get_src();
    vertex_range.second = edges.back().get_src();
  }
  return res;
}
} // namespace graphs
