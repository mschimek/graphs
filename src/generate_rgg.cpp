#include "interface.hpp"

#include "AmsSort/AmsSort.hpp"
#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "weight_generators.hpp"

namespace graphs {
std::size_t find_smaller_square_power_of_two(std::size_t i) {
  if (i <= 3)
    return 1;
  const std::size_t log_i = static_cast<std::size_t>(std::log2(i));
  const std::size_t log_i_even = log_i % 2 == 0 ? log_i : log_i - 1;
  return 1ull << log_i_even;
}

std::size_t find_next_greater_square_power_of_two(std::size_t i) {
  if (i <= 3)
    return 1;
  const std::size_t log_i = static_cast<std::size_t>(std::log2(i));
  const std::size_t log_i_even = log_i % 2 == 0 ? log_i : log_i + 1;
  return 1ull << log_i_even;
}
std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D_using_sorting(std::size_t log_n, double radius, MPIComm comm) {
  const std::size_t square_power_of_two =
      find_smaller_square_power_of_two(comm.size);
  const bool is_power_of_two_square = square_power_of_two == comm.size;
  std::pair<std::vector<WEdge>, VertexRange> result;
  if (is_power_of_two_square) {
    kagen::KaGen gen(comm.rank, comm.size);
    result = gen.Generate2DRGG(WeightGenerator<VId, Weight, WEdge>{},
                               1ull << log_n, radius);
  } else {
    if (comm.rank < square_power_of_two) {
      kagen::KaGen gen(comm.rank, square_power_of_two);
      result = gen.Generate2DRGG(WeightGenerator<VId, Weight, WEdge>{},
                                 1ull << log_n, radius);
    }
  }
  remove_upside_down(result.first, result.second);

  MPI_Datatype mpi_edge_type;
  MPI_Type_contiguous(sizeof(WEdge), MPI_BYTE, &mpi_edge_type);
  MPI_Type_commit(&mpi_edge_type);
  int tag = 100000;
  std::mt19937_64 gen(comm.rank * 100);
  auto SrcDstSort = [](const WEdge& lhs, const WEdge& rhs) {
    return std::tie(lhs.src, lhs.dst) < std::tie(rhs.src, rhs.dst);
  };
  const int num_levels = 2;
  Ams::sortLevel(mpi_edge_type, result.first, num_levels, gen, comm.comm,
                 SrcDstSort);
  MPI_Type_free(&mpi_edge_type);
  result.second.first = result.first.empty() ? -1 : result.first.front().src;
  result.second.second = result.first.empty() ? -1 : result.first.back().src;

  return result;
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, double radius, MPIComm comm) {
  const std::size_t square_power_of_two =
      find_next_greater_square_power_of_two(comm.size);
  kagen::KaGen gen(comm.rank, comm.size);
  auto result = gen.Generate2DRGG(WeightGenerator<VId, Weight, WEdge>{},
                                  1ull << log_n, radius, square_power_of_two);
  remove_upside_down(result.first, result.second);
  repair_edges(result.first, result.second, comm);
  return result;
}
} // namespace graphs
