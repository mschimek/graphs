#include <iostream>

#include "mpi.h"

#include "interface.hpp"
#include "include/utils.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  auto [edges, range] = graphs::get_rgg2D(5, 0.1);
  graphs::execute_in_order(comm, [&]() {
    for (const auto& edge : edges) {
      std::cout << edge << std::endl;
    }
  });
  MPI_Finalize();
}
