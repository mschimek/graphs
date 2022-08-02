#pragma once

#include <algorithm>
#include <cmath>
#include <thread>

#include "kagen.h"

#include "interface.hpp"
#include "mpi.h"

namespace graphs {
void remove_upside_down(WEdgeList& edges, const VertexRange& range);
void repair_edges(WEdgeList& edges, const VertexRange& range, MPIComm comm);
std::pair<std::vector<WEdge14>, VertexRange>
add_weights(UniformRandomWeightGenerator<VId, Weight, WEdge>& w_gen,
            kagen::KaGenResult&& kagen_result, MPIComm comm);

std::size_t get_next_pow_two_with_exp_divisible_by(std::size_t i,
                                                   std::size_t divisor);

template <typename F> void execute_in_order(MPIComm comm, F&& f) {
  for (std::size_t i = 0; i < comm.size; ++i) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    MPI_Barrier(comm.comm);
    if (i == comm.rank) {
      std::cout << "On rank: " << comm.rank << std::endl;
      f();
    }
  }
}

} // namespace graphs
