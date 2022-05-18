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
  std::cout << infile << std::endl;
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
  std::cout << "rank: " << comm.rank << " file size: " << file_size << " "
            << sizeof(std::fstream::pos_type)
            << " local_chunk_size: " << local_chunk_size << std::endl;
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
    std::cout << e << std::endl;
    *(out_it) = e;
    ++out_it;
  }
  return local_chunk_size;
}

template <typename Container, typename Infiles>
Container read_weighted_binaries(const Infiles& infiles,
                                 std::size_t num_VId_bytes,
                                 MPIComm comm = MPIComm{}) {
  using EdgeType = typename Container::value_type;
  std::size_t num_own_edges = 0;
  for (const auto& infile : infiles) {
    num_own_edges += compute_number_edges(infile, num_VId_bytes, comm);
  }
  Container edges;
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
           std::make_tuple(rhs.get_src(), rhs.get_dst(), lhs.get_weight());
  };
  auto SrcDstWeightEqual = [](const EdgeType& lhs, const EdgeType& rhs) {
    return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) ==
           std::make_tuple(rhs.get_src(), rhs.get_dst(), lhs.get_weight());
  };
  const int num_levels = 2;
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstWeightOrder);
  edges.erase(std::unique(edges.begin(), edges.end(), SrcDstWeightEqual), edges.end());
  Ams::sortLevel(mpi_edge_type, edges, num_levels, gen, comm.comm, SrcDstWeightOrder);
  return edges;
}
} // namespace graphs
