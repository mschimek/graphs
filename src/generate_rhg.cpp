#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rhg(std::size_t log_n, std::size_t avg_degree, double gamma, MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  auto res = gen.GenerateRHG(WeightGenerator<VId, Weight, WEdge>{},
                             1ull << log_n, gamma, avg_degree);
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  repair_edges(res.first, res.second, comm);
  return res;
}
} // namespace graphs