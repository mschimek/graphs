#include <iostream>

#include <mpi.h>

#include "../interface.hpp"
int main() {
  std::cout << "hello" << std::endl;
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //auto [edges, range] = graphs::read_unweighted_graph("/tmp/mozilla_matthias0/road-road-usa/road-road-usa.mtx", graphs::GraphFormat::MatrixMarket);
  //std::cout << rank << " " << edges.front() << " " << edges.back() << std::endl;
  auto [edges2, range2] = graphs::get_rgg2D(7,0.1);
  //std::cout << rank << " " << edges2.front() << " " << edges2.back()  << " " << edges2.size() << std::endl;
  //for(auto& edge : edges2)
  //  std::cout << edge << std::endl;
  //auto [edges3, range3] = graphs::get_rmat_edges_evenly_distributed(graphs::RMatParams(4,5));
  //std::cout << rank << " " << edges3.front() << " " << edges3.back()  << " " << edges3.size() << std::endl;
  MPI_Finalize();
}
