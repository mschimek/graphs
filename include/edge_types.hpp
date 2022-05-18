#pragma once
#include <cstdint>
#include <mpi.h>
#include <ostream>

namespace graphs {

struct WEdge14 {
  WEdge14() = default;
  WEdge14(std::uint64_t src, std::uint64_t dst, std::uint8_t weight) {
    set_src(src);
    set_dst(dst);
    set_weight(weight);
  }
  std::uint16_t src_low;
  std::uint16_t src_mid;
  std::uint16_t src_high;
  std::uint16_t dst_low;
  std::uint16_t dst_mid;
  std::uint16_t dst_high;
  std::uint8_t weight;
  std::uint8_t get_weight() const { return weight; }
  void set_weight(uint8_t weight) { this->weight = weight; }
  std::uint64_t get_src() const {
    std::uint64_t src = src_high;
    src <<= 16;
    src |= src_mid;
    src <<= 16;
    src |= src_low;
    return src;
  }
  std::uint64_t get_dst() const {
    std::uint64_t dst = dst_high;
    dst <<= 16;
    dst |= dst_mid;
    dst <<= 16;
    dst |= dst_low;
    return dst;
  }

  void set_src(uint64_t src) {
    this->src_low = src;
    this->src_mid = src >> 16;
    this->src_high = src >> 32;
  }
  void set_dst(uint64_t dst) {
    this->dst_low = dst;
    this->dst_mid = dst >> 16;
    this->dst_high = dst >> 32;
  }
  struct MPI_Type {
    MPI_Type() {
      MPI_Type_contiguous(sizeof(WEdge14), MPI_BYTE, &mpi_datatype_);
      MPI_Type_commit(&mpi_datatype_);
    }
    ~MPI_Type() { MPI_Type_free(&mpi_datatype_); }
    MPI_Datatype get_mpi_type() const { return mpi_datatype_; }
    MPI_Datatype mpi_datatype_;
  };
  friend std::ostream& operator<<(std::ostream& out, const WEdge14& edge) {
    return out << "(" << edge.get_src() << ", " << edge.get_dst() << ", "
               << std::uint32_t(edge.get_weight()) << ")";
  }
};

struct WEdge24 {
  std::uint64_t src;
  std::uint64_t dst;
  std::uint64_t weight;
  void set_src(std::uint64_t src) { this->src = src; }
  void set_dst(std::uint64_t dst) { this->dst = dst; }
  void set_weight(std::uint64_t weight) { this->weight = weight; }
  std::uint64_t get_src() const { return src; }
  std::uint64_t get_dst() const { return dst; }
  std::uint64_t get_weight() const { return weight; }
  WEdge24() = default;
  WEdge24(std::uint64_t src, std::uint64_t dst, std::uint32_t weight)
      : src{src}, dst{dst}, weight{weight} {}
  friend std::ostream& operator<<(std::ostream& out, const WEdge24& edge) {
    return out << "(" << edge.src << ", " << edge.dst << ", " << edge.weight
               << ")";
  }
  struct MPI_Type {
    MPI_Type() {
      MPI_Type_contiguous(sizeof(WEdge24), MPI_BYTE, &mpi_datatype_);
      MPI_Type_commit(&mpi_datatype_);
    }
    ~MPI_Type() { MPI_Type_free(&mpi_datatype_); }
    MPI_Datatype get_mpi_type() const { return mpi_datatype_; }
    MPI_Datatype mpi_datatype_;
  };
};
} // namespace graphs
