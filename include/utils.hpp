#pragma once

#include <algorithm>
#include <thread>

#include "mpi.h"

#include "interface.hpp"

namespace graphs {
inline void remove_upside_down(WEdgeList& edges, const VertexRange& range) {
  auto it = std::remove_if(edges.begin(), edges.end(), [&](const WEdge& edge) {
    return edge.src < range.first || edge.src > range.second;
  });
  edges.erase(it, edges.end());
}
template <typename F>
void execute_in_order(MPIComm comm, F&& f) {
  for (std::size_t i = 0; i < comm.size; ++i) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    MPI_Barrier(comm.comm);
    if (i == comm.rank) {
      std::cout << "On rank: " << comm.rank << std::endl;
      f();
    }
  }
}

}  // namespace graphs
