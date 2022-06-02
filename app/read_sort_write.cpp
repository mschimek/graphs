#include <iostream>

#include <mpi.h>
#include "interface.hpp"
#include "include/utils.hpp"

#include "tlx/cmdline_parser.hpp"

int main(int argc, char* argv[]) {
  tlx::CmdlineParser cp;
  std::vector<std::string> infiles; 
  graphs::Configs config;
  config.num_vid_bytes_in = 4;
  config.num_vid_bytes_out = 4;
  config.num_weight_bytes_in = 1;
  config.num_weight_bytes_out = 1;

  cp.set_description("tool for adding weights and backedges to unweighted, directed graph");
  cp.set_author("Matthias Schimek");
  cp.add_string("infiles", config.infile,
                "path to weighted binary infile");
  cp.add_string("outfiles", config.outfile,
                "path to weighted binary outfile");
  cp.add_size_t("num_vid_bytes_in", config.num_weight_bytes_in,
                " number of bytes for vertex ids - 4 or 8");
  cp.add_size_t("num_weight_bytes_in", config.num_weight_bytes_in,
                " number of bytes for weights");
  cp.add_size_t("num_vid_bytes_out", config.num_weight_bytes_out,
                " number of bytes for vertex ids - 4 or 8");
  cp.add_size_t("num_weight_bytes_out", config.num_weight_bytes_out,
                " number of bytes for weights");

  cp.process(argc, argv);
  MPI_Init(nullptr, nullptr);
  graphs::MPIComm comm;
  graphs::read_sort_write(config);
  MPI_Finalize();
  //graphs::read_weighted_binary(outfile, num_vid_bytes);
}

