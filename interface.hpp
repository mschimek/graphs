#pragma once

#include <cstdint>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include <mpi.h>

namespace graphs {
using VertexRange = std::pair<std::size_t, std::size_t>;
using VId = uint64_t;
using Weight = uint32_t;

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

struct WEdge {
  VId src;
  VId dst;
  Weight weight;
  WEdge() = default;
  WEdge(std::uint64_t src, std::uint64_t dst, std::uint32_t weight)
      : src{src}, dst{dst}, weight{weight} {}
  friend std::ostream& operator<<(std::ostream& out, const WEdge& edge) {
    return out << "(" << edge.src << ", " << edge.dst << ", " << edge.weight
               << ")";
  }
  struct MPI_Type {
    MPI_Type() {
      MPI_Type_contiguous(sizeof(WEdge), MPI_BYTE, &mpi_datatype_);
      MPI_Type_commit(&mpi_datatype_);
    }
    ~MPI_Type() { MPI_Type_free(&mpi_datatype_); }
    MPI_Datatype get_mpi_type() const { return mpi_datatype_; }
    MPI_Datatype mpi_datatype_;
  };
};

struct SrcDstOrder {
  bool operator()(const WEdge& lhs, const WEdge& rhs) const {
    return std::tie(lhs.src, lhs.dst) < std::tie(rhs.src, rhs.dst);
  }
};

struct SrcDstWeightOrder {
  bool operator()(const WEdge& lhs, const WEdge& rhs) const {
    return std::tie(lhs.src, lhs.dst, lhs.weight) <
           std::tie(rhs.src, rhs.dst, rhs.weight);
  }
};

inline bool operator==(const WEdge& lhs, const WEdge& rhs) {
  return std::tie(lhs.src, lhs.dst, lhs.weight) ==
         std::tie(rhs.src, rhs.dst, rhs.weight);
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

using WEdgeList = std::vector<WEdge>;

std::pair<std::vector<WEdge>, VertexRange>
get_gnm(std::size_t log_n, std::size_t log_m, MPIComm comm = MPIComm{});
std::pair<std::vector<WEdge>, VertexRange>
get_rgg2D(std::size_t log_n, double radius, MPIComm comm = MPIComm{});
std::pair<std::vector<WEdge>, VertexRange> get_rhg(std::size_t log_n,
                                                   std::size_t avg_degree,
                                                   double gamma,
                                                   MPIComm comm = MPIComm{});

std::pair<WEdgeList, VertexRange>
get_rmat_edges_evenly_distributed(const RMatParams& params,
                                  bool remove_duplicates = true,
                                  MPIComm comm = MPIComm{});

enum class GraphFormat { MatrixMarket, Snap };
std::pair<std::vector<WEdge>, VertexRange>
read_weighted_graph(const std::string& filename, GraphFormat format,
                    MPIComm comm = MPIComm{});
std::pair<std::vector<WEdge>, VertexRange>
read_unweighted_graph(const std::string& filename, GraphFormat format,
                      MPIComm comm = MPIComm{});
}; // namespace graphs



