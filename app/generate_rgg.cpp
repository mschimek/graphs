#include <iostream>
#include <algorithm>

#include "mpi.h"

#include "interface.hpp"
#include "include/utils.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  auto [edges, range] = graphs::get_rgg2D(10, 0.1);
  std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs) { return lhs.weight < rhs.weight; });
  graphs::execute_in_order(comm, [&]() {
      std::size_t diff = 0;
      std::size_t i = 1;
    for (const auto& edge : edges) {
      const auto dif = std::labs(static_cast<int>(edge.dst) - static_cast<int>(edge.src));
      diff += dif;
      std::cout << edge << " " << dif << " avg: " << (static_cast<double>(diff) / i) <<  std::endl;
      ++i;
    }
  });
  MPI_Finalize();
}
