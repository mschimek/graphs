#include <algorithm>
#include <iostream>

#include "mpi.h"

#include "include/utils.hpp"
#include "interface.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  std::size_t log_n = 15ul + std::log2(comm.size);
  std::size_t log_m = 17ul + std::log2(comm.size);
  auto [edges, range] = graphs::get_rgg3D(log_n, log_m);
  // std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs)
  // {
  //   return lhs.weight < rhs.weight;
  // });
  //
  std::sort(edges.begin(), edges.end(),
            graphs::SrcDstWeightOrder<graphs::WEdge>{});
  auto is_local = [&](const auto& v) {
    return range.first <= v && v <= range.second;
  };

  const auto num_local_edges =
      std::count_if(edges.begin(), edges.end(), [&](const auto& edge) {
        return is_local(edge.get_src()) && is_local(edge.get_dst());
      });

  graphs::execute_in_order(comm, [&]() {
    std::cout << range.first << " - " << range.second << std::endl;
    std::cout << "ratio local: "
              << static_cast<double>(num_local_edges) / edges.size()
              << std::endl;
    std::cout << edges.size() << std::endl;
    for (const auto& edge : edges) {
      std::cout << edge << std::endl;
    }
  });
  MPI_Finalize();
}
