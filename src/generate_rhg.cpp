#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rhg(std::size_t log_n, std::size_t avg_degree, double gamma,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  auto res = gen.GenerateRHG(
      UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
      1ull << log_n, gamma, avg_degree);
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder<WEdge>{});
  repair_edges(res.first, res.second, comm);
  return res;
}

std::pair<std::vector<WEdge>, VertexRange>
get_rhg_explicit_num_edges(std::size_t log_n, std::size_t log_m, double gamma,
                           WeightGeneratorConfig<Weight> wgen_config,
                           MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  const std::size_t avg_degree = (1ull << log_m) / (1ull << log_n);
  return get_rhg(log_n, avg_degree, gamma, wgen_config, comm);
}
} // namespace graphs
