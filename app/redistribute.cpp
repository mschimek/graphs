#include <iostream>
#include <random>

#include "mpi.h"

#include "include/utils.hpp"
#include "interface.hpp"

std::vector<std::size_t> choose_sizes(std::size_t m, std::size_t p) {
  if (p == 1) {
    return {m};
  }
  std::mt19937 gen(444);
  std::size_t remaining_elements = m;
  std::vector<std::size_t> sizes(p, 0);
  std::uniform_int_distribution<> distrib(0, m - 1);
  for (std::size_t i = 0; i + 1< p && remaining_elements != 0; ++i) {
    sizes[i] = distrib(gen) % remaining_elements;
    remaining_elements -= sizes[i];
  }
  sizes.back() = remaining_elements;
  std::shuffle(sizes.begin(), sizes.end(), gen);
  return sizes;
}
int main() {
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  const std::size_t log_num_edges = 18;
  auto [edges, range] = graphs::get_gnm(14ul, log_num_edges);
  const std::uint64_t num_generated_edges = edges.size();
  std::uint64_t global_m = num_generated_edges;
  MPI_Allreduce(MPI_IN_PLACE, &global_m, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  //std::cout << "global m : " << global_m << std::endl;
  const auto sizes = choose_sizes(global_m, comm.size);
  std::size_t redistribute_num = sizes[comm.rank];
  auto redistributed_edges =
      graphs::redistribute(edges, redistribute_num, comm);
  auto redistributed_edges_ =
      graphs::redistribute(redistributed_edges, num_generated_edges, comm);
  if (redistributed_edges_ != edges) {
    std::cout << "wrong res: " << comm.rank << std::endl;
    std::abort();
  }
  MPI_Finalize();
}
