#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

inline std::size_t compute_next_number_divisible_by(const std::size_t i,
                                                    const std::size_t divisor) {
  const std::size_t remainder = i % divisor;
  const std::size_t offset = remainder == 0 ? 0 : (divisor - remainder);
  return i + offset;
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid2D(std::size_t log_x, std::size_t log_y, double p, bool is_periodic,
           MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  const std::size_t square_power_of_two =
      get_next_pow_two_with_exp_divisible_by(comm.size, 2);
  std::cout << square_power_of_two << std::endl;
  auto res =
      gen.Generate2DGrid(WeightGenerator<VId, Weight, WEdge>{}, 1ull << log_x,
                         1ull << log_y, p, is_periodic, square_power_of_two);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  // graphs::execute_in_order(comm, [&]() {
  //   std::cout << res.second.first << " " << res.second.second << std::endl;
  //   for (const auto& edge : res.first) {
  //     std::cout << edge << std::endl;
  //   }
  // });
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  repair_edges(res.first, res.second, comm);
  return res;
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid2D(std::size_t log_n, double p, bool is_periodic, MPIComm comm) {
  const std::size_t log_n_even = compute_next_number_divisible_by(log_n, 2);
  return get_grid2D(log_n_even / 2, log_n_even / 2, p, is_periodic, comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid3D(std::size_t log_x, std::size_t log_y, std::size_t log_z, double p,
           bool is_periodic, MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  const std::size_t cube_power_of_two =
      get_next_pow_two_with_exp_divisible_by(comm.size, 3);
  std::cout << cube_power_of_two << std::endl;
  auto res = gen.Generate3DGrid(WeightGenerator<VId, Weight, WEdge>{},
                                1ull << log_x, 1ull << log_y, 1ull << log_z, p,
                                is_periodic, cube_power_of_two);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  //graphs::execute_in_order(comm, [&]() {
  //  std::cout << res.second.first << " " << res.second.second << std::endl;
  //  for (const auto& edge : res.first) {
  //    std::cout << edge << std::endl;
  //  }
  //});
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  repair_edges(res.first, res.second, comm);
  return res;
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid3D(std::size_t log_n, double p, bool is_periodic, MPIComm comm) {
  const std::size_t log_n_div_by_three =
      compute_next_number_divisible_by(log_n, 3);
  return get_grid3D(log_n_div_by_three / 3, log_n_div_by_three / 3,
                    log_n_div_by_three / 3, p, is_periodic, comm);
}
} // namespace graphs
