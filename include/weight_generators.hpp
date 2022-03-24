#pragma once

#include <algorithm>
#include <cstdint>
#include <type_traits>

namespace graphs {
template <typename Vertex, typename Weight, typename Edge>
struct WeightGenerator {
  using EdgeType = Edge;
  using WeightType = Weight;
  WeightGenerator(Weight min_weight, Weight max_weight) : min_weight_{min_weight}, max_weight_{max_weight} {}
  WeightGenerator() : WeightGenerator(1, 100'000) {}
  Weight operator()(Vertex src, Vertex dst) {
    return get_random_weight(src, dst, min_weight_, max_weight_);
  }
  Weight get_min_weight() const { return min_weight_; }
  Weight get_max_weight() const { return max_weight_; }

 private:
  Weight min_weight_;
  Weight max_weight_;

  // from https://stackoverflow.com/a/12996028
  std::uint32_t hash32(std::uint32_t x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }
  std::uint64_t hash64(std::uint64_t x) {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
  }
  std::size_t knuth(std::size_t x) { return x * 2654435761 % (1ull << 31); }

  Weight get_random_weight(const Vertex& src, const Vertex& dst,
                           const Weight min_weight = 1,
                           const Weight max_weight = 100'000) {
    static_assert(std::is_integral_v<Vertex>,
                  "Vertex type must be an integer.");
    static_assert(std::is_integral_v<Weight>,
                  "Weight type must be an integer.");
    Weight res = 0;
    if constexpr (sizeof(Vertex) <= 4) {
      unsigned int lower_part_src = static_cast<unsigned int>(src);
      unsigned int lower_part_dst = static_cast<unsigned int>(dst);
      res = hash32((hash32(src + 111) + hash32(dst + 111)));
    } else {
      unsigned int upper_part_src = static_cast<unsigned int>(src >> 32);
      unsigned int upper_part_dst = static_cast<unsigned int>(dst >> 32);
      res = hash64(hash64(src + 111) + hash64(dst + 111));
    }
    return std::max((res % max_weight), min_weight);
  }
};
}  // namespace graphs
