#include <iostream>

#include <mpi.h>
#include "interface.hpp"
#include "include/utils.hpp"

#include "tlx/cmdline_parser.hpp"

int main(int argc, char* argv[]) {
  tlx::CmdlineParser cp;
  std::vector<std::string> infiles; 
  std::size_t num_vid_bytes = 8;
  std::size_t num_weight_bytes = 1;
  cp.set_description("tool for adding weights and backedges to unweighted, directed graph");
  cp.set_author("Matthias Schimek");
  cp.add_stringlist("infiles", infiles,
                "path to unweighted, directed graph");
  cp.add_size_t("num_vid_bytes", num_vid_bytes,
                " number of bytes for vertex ids - 4 or 8");
  cp.add_size_t("num_weight_bytes", num_weight_bytes,
                " number of bytes for weights");

  cp.process(argc, argv);
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  num_vid_bytes = num_vid_bytes == 4 ? 4 : 8;
  auto edges = graphs::read_weighted_binaries<graphs::WEdgeList, std::vector<std::string>>(infiles, num_vid_bytes, num_weight_bytes, comm);
  graphs::execute_in_order(comm, [&](){ for(const auto& elem : edges.first) { std::cout << elem << std::endl;}});
  MPI_Finalize();
  //graphs::read_weighted_binary(outfile, num_vid_bytes);
}

