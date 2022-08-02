#include <iostream>

#include "mpi.h"

#include "include/utils.hpp"
#include "interface.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  auto [edges, range] = graphs::get_gnm(8ul, 10ul);
  graphs::execute_in_order(comm, [&]() {
    std::cout << "num edges: " << edges.size() << std::endl;
     for (const auto& edge : edges) {
       std::cout << edge << std::endl;
     }
  });
  MPI_Finalize();
}
