#include "interface.hpp"

#include <algorithm>

#include "kagen.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"
#include "memory_utils.hpp"

namespace graphs {

std::pair<std::vector<WEdge14>, VertexRange>
get_gnm(std::size_t log_n, std::size_t log_m,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  memory_stats().print("before kagen 1");
  kagen::KaGen gen(MPI_COMM_WORLD);
  memory_stats().print("before kagen 2");
  auto res = gen.GenerateUndirectedGNM(1ull << log_n, 1ull << log_m);
  UniformRandomWeightGenerator<VId, Weight, WEdge> w_gen{wgen_config};
  return add_weights(w_gen, std::move(res), comm);
}
} // namespace graphs
