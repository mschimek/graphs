#include <algorithm>
#include <iostream>

#include "mpi.h"

#include "include/utils.hpp"
#include "interface.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  auto [edges, range] = graphs::get_rhg(12, 3, 8);
  // std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs)
  // {
  //   return lhs.weight < rhs.weight;
  // });
  //
  std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs) {
    return std::tie(lhs.src, lhs.dst) < std::tie(rhs.src, rhs.dst);
  });
  auto is_local = [&](const auto& v) {
    return range.first <= v && v <= range.second;
  };

  const auto num_local_edges =
      std::count_if(edges.begin(), edges.end(), [&](const auto& edge) {
        return is_local(edge.src) && is_local(edge.dst);
      });

  graphs::execute_in_order(comm, [&]() {
    std::cout << range.first << " - " << range.second << std::endl;
    std::cout << "ratio local: "
              << static_cast<double>(num_local_edges) / edges.size()
              << std::endl;
    std::cout << edges.size() << std::endl;
  });
  MPI_Finalize();
}
