#pragma once

#include <algorithm>
#include <thread>

#include "interface.hpp"
#include "mpi.h"

namespace graphs {
void remove_upside_down(WEdgeList& edges, const VertexRange& range);
void repair_edges(WEdgeList& edges, const VertexRange& range, MPIComm comm);

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
