#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::size_t find_next_greater_square_power_of_two(std::size_t i) {
  if (i <= 3)
    return 1;
  const std::size_t log_i = static_cast<std::size_t>(std::log2(i));
  const std::size_t log_i_even = log_i % 2 == 0 ? log_i : log_i + 1;
  return 1ull << log_i_even;
}

std::pair<std::vector<WEdge>, VertexRange> get_grid(std::size_t log_x,
                                                    std::size_t log_y, double p,
                                                    bool is_periodic,
                                                    MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  const std::size_t square_power_of_two =
      find_next_greater_square_power_of_two(comm.size);
  auto res =
      gen.Generate2DGrid(WeightGenerator<VId, Weight, WEdge>{}, 1ull << log_x,
                         1ull << log_y, p, is_periodic, square_power_of_two);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  graphs::execute_in_order(comm, [&]() {
    std::cout << res.second.first << " " << res.second.second << std::endl;
    for (const auto& edge : res.first) {
      std::cout << edge << std::endl;
    }
  });
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  repair_edges(res.first, res.second, comm);
  return res;
}
} // namespace graphs
