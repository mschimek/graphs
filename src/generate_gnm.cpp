#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"
#include "ips4o/ips4o.hpp"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"
#include "memory_utils.hpp"

namespace graphs {

std::pair<std::vector<WEdge14>, VertexRange>
get_gnm(std::size_t log_n, std::size_t log_m,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  memory_stats().print("before kagen 1");
  kagen::KaGen gen(comm.rank, comm.size);
  memory_stats().print("before kagen 2");
  auto res = gen.GenerateUndirectedGNM(
      UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
      1ull << log_n, 1ull << log_m);
  memory_stats().print("after getting  edges");
  remove_upside_down(res.first, res.second);
  memory_stats().print("after remove upside down");
  ips4o::parallel::sort(res.first.begin(), res.first.end(), SrcDstOrder<WEdge14>{});
  repair_edges(res.first, res.second, comm);
  res.first.shrink_to_fit();
  memory_stats().print("after repair");
  return res;
}
} // namespace graphs
