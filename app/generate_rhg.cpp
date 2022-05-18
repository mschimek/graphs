#include <algorithm>
#include <iostream>

#include "mpi.h"

#include "include/utils.hpp"
#include "interface.hpp"

int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;

  std::uint64_t logn = 0;
  std::uint64_t avg_deg = 0;
  double power_law = 0;
  if (comm.rank == 0) {
    std::cout << "enter// logn avg_deg powerlaw_exp" << std::endl;
    std::cin >> logn >> avg_deg >> power_law;
  }
  MPI_Allreduce(MPI_IN_PLACE, &logn, 1, MPI_UINT64_T, MPI_MAX, comm.comm);
  MPI_Allreduce(MPI_IN_PLACE, &avg_deg, 1, MPI_UINT64_T, MPI_MAX, comm.comm);
  MPI_Allreduce(MPI_IN_PLACE, &power_law, 1, MPI_DOUBLE, MPI_MAX, comm.comm);

  auto [edges, range] = graphs::get_rhg(logn, avg_deg, power_law);
  std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs) {
    return std::make_pair(lhs.get_src(), lhs.get_dst()) < std::make_pair(rhs.get_src(), rhs.get_dst());
  });
  auto is_local = [&](const auto& v) {
    return range.first <= v && v <= range.second;
  };

  const auto num_local_edges =
      std::count_if(edges.begin(), edges.end(), [&](const auto& edge) {
        return is_local(edge.get_src()) && is_local(edge.get_dst());
      });

  constexpr bool output_edges = false;

  graphs::execute_in_order(comm, [&]() {
    std::cout << "vertex range: " << range.first << " - " << range.second
              << std::endl;
    std::cout << "ratio local: "
              << static_cast<double>(num_local_edges) / edges.size()
              << std::endl;
    std::cout << "num edges: " << edges.size() << std::endl;
    if (output_edges) {
      for (const auto& edge : edges) {
        std::cout << edge << std::endl;
      }
    }
  });
  MPI_Finalize();
}
