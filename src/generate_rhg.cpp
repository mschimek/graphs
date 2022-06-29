#include "interface.hpp"

#include <algorithm>

#include "KaGen/interface/kagen_interface.h"

#include "include/utils.hpp"
#include "include/weight_generators.hpp"

namespace graphs {

std::pair<std::vector<WEdge>, VertexRange>
get_rhg(std::size_t log_n, std::size_t avg_degree, double gamma,
        WeightGeneratorConfig<Weight> wgen_config, MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  auto res = gen.GenerateRHG(
      UniformRandomWeightGenerator<VId, Weight, WEdge>{wgen_config},
      1ull << log_n, gamma, avg_degree);
  std::uint64_t num_edges = res.first.size();
  MPI_Allreduce(MPI_IN_PLACE, &num_edges, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  if (comm.rank == 0) {
    std::cout << "num edges before remove upside down: " << num_edges
              << std::endl;
  }
  remove_upside_down(res.first, res.second);
  num_edges = res.first.size();
  MPI_Allreduce(MPI_IN_PLACE, &num_edges, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  if (comm.rank == 0) {
    std::cout << "num edges after remove upside down: " << num_edges
              << std::endl;
  }
  std::sort(res.first.begin(), res.first.end(), SrcDstOrder<WEdge>{});
  repair_edges(res.first, res.second, comm);
  num_edges = res.first.size();
  MPI_Allreduce(MPI_IN_PLACE, &num_edges, 1, MPI_UINT64_T, MPI_SUM, comm.comm);
  if (comm.rank == 0) {
    std::cout << "num edges after repair: " << num_edges << std::endl;
  }
  return res;
}

std::pair<std::vector<WEdge>, VertexRange>
get_rhg_explicit_num_edges(std::size_t log_n, std::size_t log_m, double gamma,
                           WeightGeneratorConfig<Weight> wgen_config,
                           MPIComm comm) {
  kagen::KaGen gen(comm.rank, comm.size);
  const std::size_t avg_degree = 2ull * (1ull << log_m) / (1ull << log_n);
  return get_rhg(log_n, avg_degree, gamma, wgen_config, comm);
}
} // namespace graphs
