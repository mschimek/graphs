#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_gnm(std::size_t log_n, std::size_t log_m,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  auto res = gen.GenerateUndirectedGNM(
      UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
      1ull << log_n, 1ull << log_m);
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  repair_edges(res.first, res.second, comm);
  return res;
}
} // namespace graphs
