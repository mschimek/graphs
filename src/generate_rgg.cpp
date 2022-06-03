#include <math.h>

#include "interface.hpp"

#include "AmsSort/AmsSort.hpp"
#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, double radius,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const std::size_t square_power_of_two =
      get_next_pow_two_with_exp_divisible_by(comm.size, 2);
  kagen::KaGen gen(comm.rank, comm.size);
  auto result = gen.Generate2DRGG(UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
                                  1ull << log_n, radius, square_power_of_two);
  remove_upside_down(result.first, result.second);
  repair_edges(result.first, result.second, comm);
  return result;
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg3D(std::size_t log_n, double radius,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const std::size_t cubic_power_of_two =
      get_next_pow_two_with_exp_divisible_by(comm.size, 3);
  kagen::KaGen gen(comm.rank, comm.size);
  auto result = gen.Generate3DRGG(UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
                                  1ull << log_n, radius, cubic_power_of_two);
  remove_upside_down(result.first, result.second);
  repair_edges(result.first, result.second, comm);
  return result;
}

double compute_radius_2d(std::size_t log_n, std::size_t log_m) {
  const double n = 1ull << log_n;
  const double m = 1ull << log_m;
  const double radius = (1 / n) * std::sqrt(2*m / M_PI);
  return radius;
}

double compute_radius_3d(std::size_t log_n, std::size_t log_m) {
  const double n = 1ull << log_n;
  const double m = 1ull << log_m;
  const double radiant = 3.0/2.0 * m / (M_PI * n * n);
  const double radius = std::cbrt(radiant);
  return radius;
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, std::size_t log_m,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const double radius = compute_radius_2d(log_n, log_m);
  return get_rgg2D(log_n, radius, wgen_config, comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg3D(std::size_t log_n, std::size_t log_m,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  const double radius = compute_radius_3d(log_n, log_m);
  return get_rgg3D(log_n, radius, wgen_config, comm);
}
} // namespace graphs
