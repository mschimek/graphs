#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <type_traits>

namespace graphs {

template <typename Weight> struct WeightGeneratorConfig {
  WeightGeneratorConfig() = default;
  Weight min_weight = 1;
  Weight max_weight = 100'000;
  bool use_random_weights = true;
  std::size_t seed = 1;
};

template <typename Vertex, typename Weight, typename Edge>
struct WeightGenerator {
  using EdgeType = Edge;
  using WeightType = Weight;
  WeightGenerator(Weight min_weight, Weight max_weight)
      : config_{min_weight, max_weight, true} {}
  WeightGenerator(Weight max_weight) : config_{1, max_weight, true} {}
  WeightGenerator(WeightGeneratorConfig<Weight> config) : config_{config} {}
  WeightGenerator() : WeightGenerator(1, 100'000) {}

  Weight operator()(Vertex src, Vertex dst) const {
    return get_random_weight(src, dst, config_.min_weight, config_.max_weight);
  }

  Weight operator()(Vertex src, Vertex dst, double dist_factor) const {
    if (config_.use_random_weights) {
      return get_random_weight(src, dst, config_.min_weight,
                               config_.max_weight);
    }
    const Weight w = std::min(
        static_cast<Weight>(dist_factor * get_max_weight()), get_max_weight());
    return std::max(get_min_weight(), w);
  }
  Weight get_min_weight() const { return config_.min_weight; }
  Weight get_max_weight() const { return config_.max_weight; }

private:
  WeightGeneratorConfig<Weight> config_;

  // from https://stackoverflow.com/a/12996028
  std::uint32_t hash32(std::uint32_t x) const {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }
  std::uint64_t hash64(std::uint64_t x) const {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
  }
  std::size_t knuth(std::size_t x) { return x * 2654435761 % (1ull << 31); }

  Weight get_random_weight(const Vertex& src, const Vertex& dst,
                           const Weight min_weight,
                           const Weight max_weight) const {
    static_assert(std::is_integral_v<Vertex>,
                  "Vertex type must be an integer.");
    static_assert(std::is_integral_v<Weight>,
                  "Weight type must be an integer.");
    Weight res = 0;
    const Vertex offset = 111;
    if constexpr (sizeof(Vertex) <= 4) {
      res = hash32((hash32(src + offset) + hash32(dst + offset)));
    } else {
      res = hash64(hash64(src + offset) + hash64(dst + offset));
    }
    return std::max((res % max_weight), min_weight);
  }
};

template <typename Vertex, typename Weight, typename Edge>
struct UniformRandomWeightGenerator {
  using EdgeType = Edge;
  using WeightType = Weight;
  UniformRandomWeightGenerator(WeightGeneratorConfig<Weight> config) : config_{config}, gen{config_.seed}, distribution(config_.min_weight, config_.max_weight) {}

  Weight operator()(Vertex src, Vertex dst) {
      return distribution(gen);
  }

  Weight operator()(Vertex src, Vertex dst, double dist_factor) {
    if (config_.use_random_weights) {
      return distribution(gen);
    }
    const Weight w = std::min(
        static_cast<Weight>(dist_factor * get_max_weight()), get_max_weight());
    return std::max(get_min_weight(), w);
  }
  Weight get_min_weight() const { return config_.min_weight; }
  Weight get_max_weight() const { return config_.max_weight; }

private:
  WeightGeneratorConfig<Weight> config_;
  std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<Weight> distribution;
};
} // namespace graphs
