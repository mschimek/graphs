#pragma once

#include "include/edge_types.hpp"
#include "include/weight_generators.hpp"
#include <cstdint>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include <mpi.h>

namespace graphs {
using VertexRange = std::pair<std::size_t, std::size_t>;
using VId = std::uint64_t;
using Weight = std::uint8_t;

struct MPIComm {
  MPIComm(MPI_Comm comm_) : comm{comm_} {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  }
  MPIComm() : MPIComm(MPI_COMM_WORLD) {}
  MPI_Comm comm;
  int rank;
  int size;
};

using WEdge = WEdge14;
using WEdgeList = std::vector<WEdge>;

template <typename EdgeType> struct SrcDstOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::make_pair(lhs.get_src(), lhs.get_dst()) <
           std::make_pair(rhs.get_src(), rhs.get_dst());
  }
};

template <typename EdgeType> struct DstSrcOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::make_tuple(lhs.get_dst(), lhs.get_src()) <
           std::make_tuple(rhs.get_dst(), rhs.get_src());
  }
};

template <typename EdgeType> struct DstSrcWeightOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::make_tuple(lhs.get_dst(), lhs.get_src(), lhs.get_weight()) <
           std::make_tuple(rhs.get_dst(), rhs.get_src(), rhs.get_weight());
  }
};

template <typename EdgeType> struct SrcDstWeightOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) <
           std::make_tuple(rhs.get_src(), rhs.get_dst(), rhs.get_weight());
  }
};

template <typename EdgeType> struct SrcDstWeightEqual {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) ==
           std::make_tuple(rhs.get_src(), rhs.get_dst(), rhs.get_weight());
  }
};

inline bool operator==(const WEdge& lhs, const WEdge& rhs) {
  return std::make_tuple(lhs.get_src(), lhs.get_dst(), lhs.get_weight()) ==
         std::make_tuple(rhs.get_src(), rhs.get_dst(), rhs.get_weight());
}

struct RMatParams {
  RMatParams(std::size_t logn, std::size_t logm) : log_n(logn), log_m(logm) {}
  RMatParams(std::size_t logn, std::size_t logm, double a, double b, double c)
      : log_n(logn), log_m(logm), a(a), b(b), c(c) {}
  std::size_t log_n = 1;
  std::size_t log_m = 1;
  std::size_t depth = log_n - 1;
  double a = 0.57, b = 0.19, c = 0.19;
};


std::pair<std::vector<WEdge>, VertexRange> get_gnm(
    std::size_t log_n, std::size_t log_m,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_rgg2D(
    std::size_t log_n, double radius,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_rgg2D(
    std::size_t log_n, std::size_t log_m,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_rhg(
    std::size_t log_n, std::size_t avg_degree, double gamma,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_rhg_explicit_num_edges(
    std::size_t log_n, std::size_t log_m, double gamma,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_grid2D(
    std::size_t log_x, std::size_t log_y, double p, bool is_periodic,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_grid2D(
    std::size_t log_n, double p, bool is_periodic,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});
std::pair<std::vector<WEdge>, VertexRange> get_grid3D(
    std::size_t log_x, std::size_t log_y, std::size_t log_z, double p,
    bool is_periodic,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<std::vector<WEdge>, VertexRange> get_grid3D(
    std::size_t log_n, double p, bool is_periodic,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    MPIComm comm = MPIComm{});

std::pair<WEdgeList, VertexRange> get_rmat_edges_evenly_distributed(
    const RMatParams& params,
    WeightGeneratorConfig<Weight> wgen_config = WeightGeneratorConfig<Weight>{},
    bool remove_duplicates = true, MPIComm comm = MPIComm{});

enum class GraphFormat { MatrixMarket, Snap, Binary };
std::pair<std::vector<WEdge>, VertexRange>
read_weighted_graph(const std::string& filename, GraphFormat format,
                    MPIComm comm = MPIComm{});
std::pair<std::vector<WEdge>, VertexRange>
read_unweighted_graph(const std::string& filename, GraphFormat format,
                      MPIComm comm = MPIComm{});

void add_weight_and_back_edges_in_directed_graph(const std::string& infile,
                               const std::string& outfile, Weight max_weight,
                               std::size_t num_VId_bytes);
void add_weight_and_back_edges_in_undirected_graph(const std::string& infile,
                               const std::string& outfile, Weight max_weight,
                               std::size_t num_VId_bytes);
}; // namespace graphs

#include "include/template_functions.hpp"
