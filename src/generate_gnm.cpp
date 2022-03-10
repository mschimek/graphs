#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/weight_generators.hpp"

namespace graphs {

void remove_upside_down(WEdgeList& edges, const VertexRange& range) {
  auto it = std::remove_if(edges.begin(), edges.end(), [&](const WEdge& edge) {
    return edge.src < range.first || edge.src > range.second;
  });
  edges.erase(it, edges.end());
}
std::pair<std::vector<WEdge>, VertexRange> get_gnm(std::size_t log_n,
                                                   std::size_t log_m,
                                                   MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  auto res = gen.GenerateUndirectedGNM(WeightGenerator<VId, Weight, WEdge>{}, 1ull << log_n,
                                       1ull << log_m);
  remove_upside_down(res.first, res.second);
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder{});
  return res;
}
}
