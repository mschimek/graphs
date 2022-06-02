#include <iostream>

#include "interface.hpp"
#include <mpi.h>

#include "tlx/cmdline_parser.hpp"

int main(int argc, char* argv[]) {
  tlx::CmdlineParser cp;
  std::string infile = "in";
  std::string outfile = "out";
  std::size_t max_weight = 254;
  std::size_t num_vid_bytes = 8;
  std::size_t num_weight_bytes = 1;
  bool is_graph_directed = false;
  bool is_graph_weighted = false;
  bool sort_output = false;

  cp.set_description(
      "tool for adding weights and backedges to unweighted, directed graph");
  cp.set_author("Matthias Schimek");
  cp.add_string("infile", infile, "path to unweighted, directed graph");
  cp.add_string("outfile", outfile, "path to outfile");
  cp.add_size_t("max_weight", max_weight,
                " max weight - needs to be smaller than 256");
  cp.add_size_t("num_vid_bytes", num_vid_bytes,
                " number of bytes for vertex ids - 4 or 8");
  cp.add_size_t("num_weight_bytes", num_weight_bytes,
                " number of bytes for weight - usually 1");
  cp.add_bool(
      "graph_is_directed", is_graph_directed,
      " directed graph == there is either edge (u,v) or (v,u) in the graph");
  cp.add_bool("graph_is_weighted", is_graph_weighted, " no weights are added");
  cp.add_bool("sort_output", sort_output, " sort edges in memory (only applicable for small files)");

  cp.process(argc, argv);
  max_weight = std::min(max_weight, std::size_t(255));
  num_vid_bytes = num_vid_bytes == 4 ? 4 : 8;
  if (is_graph_weighted) {
    graphs::transform_undirected_graph_to_binary(
        infile, outfile, max_weight, num_vid_bytes, num_weight_bytes, sort_output);
  } else {
    if (is_graph_directed) {
      graphs::add_weight_and_back_edges_in_directed_graph(
          infile, outfile, max_weight, num_vid_bytes, num_weight_bytes);
    } else {
      graphs::add_weight_and_back_edges_in_undirected_graph(
          infile, outfile, max_weight, num_vid_bytes, num_weight_bytes);
    }
  }
  // graphs::read_weighted_binary(outfile, num_vid_bytes);
}
