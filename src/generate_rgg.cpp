#include <math.h>

#include "interface.hpp"

#include "kagen.h"

#include "include/utils.hpp"
#include "weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, double radius,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRGG2D(1ull << log_n, radius);
  return add_weights(wgen, std::move(result), comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg3D(std::size_t log_n, double radius,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRGG3D(1ull << log_n, radius);
  return add_weights(wgen, std::move(result), comm);

}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, std::size_t log_m,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRGG2D_NM(1ull << log_n, 1ull << log_m);
  return add_weights(wgen, std::move(result), comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_rgg3D(std::size_t log_n, std::size_t log_m,
          WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRGG3D_NM(1ull << log_n, 1ull << log_m);
  return add_weights(wgen, std::move(result), comm);
}
} // namespace graphs
