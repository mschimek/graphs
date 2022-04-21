#include "interface.hpp"

#include "AmsSort/AmsSort.hpp"
#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, double radius, MPIComm comm) {
  const std::size_t square_power_of_two =
      get_next_pow_two_with_exp_divisible_by(comm.size, 2);
  kagen::KaGen gen(comm.rank, comm.size);
  auto result = gen.Generate2DRGG(WeightGenerator<VId, Weight, WEdge>{},
                                  1ull << log_n, radius, square_power_of_two);
  remove_upside_down(result.first, result.second);
  repair_edges(result.first, result.second, comm);
  return result;
}
} // namespace graphs
