#pragma once

#include <cstdint>

namespace graphs {
template <typename Vertex, typename Weight, typename Edge>
struct WeightGenerator {
  using EdgeType = Edge;
  Weight operator()(Vertex src, Vertex dst) {
    return get_random_weight(src, dst);
  }

 private:
  unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }
  std::size_t knuth(std::size_t x) { return x * 2654435761 % (1ull << 31); }

  Weight get_random_weight(const Vertex& src, const Vertex& dst,
                           const Weight max_weight = 100'000) {
    Weight res = 0;
    unsigned int lower_part_src = static_cast<unsigned int>(src);
    unsigned int lower_part_dst = static_cast<unsigned int>(dst);
    res = (hash(lower_part_src) + hash(lower_part_dst));
    if constexpr (sizeof(Vertex) > 4) {
      unsigned int upper_part_src = static_cast<unsigned int>(src >> 32);
      unsigned int upper_part_dst = static_cast<unsigned int>(dst >> 32);
      res += hash(upper_part_src) + hash(upper_part_dst);
    }
    return (res % max_weight) +
           1u;  // weight 0 is not allowed for input graph edges
  }
};
}  // namespace graph
