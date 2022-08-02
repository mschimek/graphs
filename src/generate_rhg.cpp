#include "interface.hpp"

#include <algorithm>

#include "kagen.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rhg(std::size_t log_n, std::size_t avg_degree, double gamma,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRHG(gamma, 1ull << log_n, avg_degree);
  return add_weights(wgen, std::move(result), comm);
}

std::pair<std::vector<WEdge>, VertexRange>
get_rhg_explicit_num_edges(std::size_t log_n, std::size_t log_m, double gamma,
                           WeightGeneratorConfig<Weight> wgen_config,
                           MPIComm comm) {
  UniformRandomWeightGenerator<VId, Weight, WEdge> wgen{wgen_config};
  kagen::KaGen gen(comm.comm);
  auto result = gen.GenerateRHG_NM(gamma, 1ull << log_n, 1ull << log_m);
  return add_weights(wgen, std::move(result), comm);
}
} // namespace graphs
