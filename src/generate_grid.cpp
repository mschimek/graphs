#include "interface.hpp"

#include <algorithm>

#include "kagen.h"

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
           WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result =
      gen.GenerateGrid2D(1ull << log_x, 1ull << log_y, p, is_periodic);
  return add_weights(wgen, std::move(result), comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid2D(std::size_t log_n, double p, bool is_periodic,
           WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const std::size_t log_n_even = compute_next_number_divisible_by(log_n, 2);
  return get_grid2D(log_n_even / 2, log_n_even / 2, p, is_periodic, wgen_config,
                    comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid3D(std::size_t log_x, std::size_t log_y, std::size_t log_z, double p,
           bool is_periodic, WeightGeneratorConfig<Weight> wgen_config,
           MPIComm comm) {

  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateGrid3D(1ull << log_x, 1ull << log_y, 1ull << log_z,
                                   p, is_periodic);
  return add_weights(wgen, std::move(result), comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_grid3D(std::size_t log_n, double p, bool is_periodic,
           WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const std::size_t log_n_div_by_three =
      compute_next_number_divisible_by(log_n, 3);
  return get_grid3D(log_n_div_by_three / 3, log_n_div_by_three / 3,
                    log_n_div_by_three / 3, p, is_periodic, wgen_config, comm);
}
} // namespace graphs
