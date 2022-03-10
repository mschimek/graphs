#include "interface.hpp"

#include <fstream>
#include <sstream>

#include "RQuick/RQuick.hpp"

#include "../include/weight_generators.hpp"

namespace graphs::MatrixMarket {
inline std::tuple<std::size_t, std::size_t, std::size_t> read_header(
    std::ifstream& in) {
  in.seekg(0, std::ios::beg);
  std::string line;
  std::size_t header_chars = 0;
  while (std::getline(in, line, '\n')) {
    header_chars += line.size() + 1;
    if (line.front() != '%') break;
  }
  std::stringstream ss(line);
  std::size_t rows = 0, columns = 0, nonzeros = 0;
  ss >> rows >> columns >> nonzeros;

  auto res = std::make_tuple(rows, nonzeros, header_chars);
  return res;
}

inline WEdgeList read_content(std::ifstream& in, std::size_t start_pos,
                              std::size_t end_pos, MPIComm comm) {
  static WeightGenerator<VId, Weight, WEdge> weight_gen;
  std::string line;
  std::size_t num_chars_read = 0;
  WEdgeList edges;
  in.seekg(start_pos - 4);
  std::getline(in, line, '\n');
  const bool do_not_discard_first_line =
      (in.tellg() == start_pos) || (comm.rank == 0);
  in.seekg(start_pos);

  if (!do_not_discard_first_line) {
    std::getline(in, line, '\n');
  }
  while (std::getline(in, line, '\n')) {
    std::stringstream ss(line);
    VId src, dst;
    ss >> src >> dst;
    const Weight w = weight_gen(src, dst);
    edges.emplace_back(src, dst, w);
    edges.emplace_back(dst, src, w);
    num_chars_read += line.size() + 1;
    if (in.tellg() >= end_pos) break;
  }
  return edges;
}

inline std::pair<WEdgeList, VertexRange> read_from_file(
    const std::string& filename, MPIComm comm) {
  std::ifstream in{filename.c_str(), std::ios::binary | std::ios::in};

  in.seekg(0, std::ios::end);
  const std::size_t length = in.tellg();
  const auto [n, m, header_chars] = read_header(in);
  std::cout << comm.rank << " length: " << length
            << " header_chars: " << header_chars << std::endl;
  in.seekg(header_chars);

  const bool is_last_rank = comm.rank + 1 == comm.size;
  const std::size_t size_without_header = length - header_chars;
  const std::size_t chunk_size = size_without_header / comm.size;
  const std::size_t chunk_remainder =
      is_last_rank ? size_without_header % comm.size : 0;
  const std::size_t start_pos = header_chars + comm.rank * chunk_size;
  const std::size_t end_pos =
      header_chars + (comm.rank + 1) * chunk_size + chunk_remainder;

  WEdgeList edges = read_content(in, start_pos, end_pos, comm);

  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(WEdge), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  int tag = 100000;
  std::mt19937_64 gen(comm.rank * 100);
  auto SrcDstSort = [](const WEdge& lhs, const WEdge& rhs) {
    return std::tie(lhs.src, lhs.dst) < std::tie(rhs.src, rhs.dst);
  };
  RQuick::sort(mpi_edge_type, edges, tag, gen, comm.comm, SrcDstSort);
  MPI_Type_free(&mpi_edge_type);
  VId v_min = edges.empty() ? -1 : edges.front().src;
  VId v_max = edges.empty() ? -1 : edges.back().src;
  return std::make_pair(std::move(edges), VertexRange{v_min, v_max});
}
}  // namespace graphs::MatrixMarket

namespace graphs {
std::pair<std::vector<WEdge>, VertexRange> read_unweighted_graph(
    const std::string& filename, GraphFormat format, MPIComm comm) {
  switch (format) {
    case GraphFormat::MatrixMarket:
      return MatrixMarket::read_from_file(filename, comm);
    default:
      return std::make_pair(WEdgeList(1), VertexRange{-1, -1});
  }
}
}  // namespace graphs
