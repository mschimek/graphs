#pragma once

#include <algorithm>

#include "interface.hpp"

namespace graphs {
inline void remove_upside_down(WEdgeList& edges, const VertexRange& range) {
  auto it = std::remove_if(edges.begin(), edges.end(), [&](const WEdge& edge) {
    return edge.src < range.first || edge.src > range.second;
  });
  edges.erase(it, edges.end());
}
}  // namespace graphs
