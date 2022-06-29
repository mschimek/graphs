#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <type_traits>

namespace graphs {

enum class DistanceType { Random, Euclidean, SquaredEuclidean };

template <typename Weight> struct WeightGeneratorConfig {
  WeightGeneratorConfig() = default;
  Weight min_weight = 1;
  Weight max_weight = 254;
  DistanceType distance_type = DistanceType::Random;
  std::size_t seed = 1;
};

template <typename Vertex, typename Weight, typename Edge>
struct WeightGenerator {
  using EdgeType = Edge;
  using WeightType = Weight;
  WeightGenerator(Weight min_weight, Weight max_weight)
      : config_{min_weight, max_weight, DistanceType::Random} {}
  WeightGenerator(Weight max_weight)
      : config_{1, max_weight, DistanceType::Random} {}
  WeightGenerator(WeightGeneratorConfig<Weight> config) : config_{config} {}
  WeightGenerator() : WeightGenerator(1, 254) {}

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
    std::uint64_t res = 0;
    const Vertex offset = 111;
    if constexpr (sizeof(Vertex) <= 4) {
      res = hash32((hash32(src + offset) + hash32(dst + offset)));
    } else {
      res = hash64(hash64(src + offset) + hash64(dst + offset));
    }
    res = std::max((res % max_weight), static_cast<std::uint64_t>(min_weight));
    return static_cast<Weight>(res);
  }
};

template <typename Vertex, typename Weight, typename Edge>
struct UniformRandomWeightGenerator {
  using EdgeType = Edge;
  using WeightType = Weight;
  UniformRandomWeightGenerator(WeightGeneratorConfig<Weight> config)
      : config_{config}, gen{config_.seed},
        distribution(config_.min_weight, config_.max_weight) {}

  Weight operator()(Vertex src, Vertex dst) { return distribution(gen); }

  Weight operator()(Vertex src, Vertex dst, double squared_euclidean_dist,
                    double radius) {
    switch (config_.distance_type) {
    case DistanceType::Random: {

      // if (rank == 0) {
      //   std::cout << "random: " << std::endl;
      // }
      return distribution(gen);
    }
    case DistanceType::Euclidean: {
      const double distance = std::sqrt(squared_euclidean_dist);
      const double dist_factor = distance / radius;
      const auto w = get_weight(dist_factor);
      // if (rank == 0) {
      //   std::cout << "euclidean: " << distance << " factor: " << dist_factor
      //             << " radius: " << radius << " w: " << int(w) << std::endl;
      // }
      return w;
    }
    case DistanceType::SquaredEuclidean: {
      const double dist_factor = squared_euclidean_dist / (radius * radius);
      const auto w = get_weight(dist_factor);
      // if (rank == 0) {
      //   std::cout << "sq_euclidean: " << squared_euclidean_dist
      //             << " factor: " << dist_factor << " radius: " << radius
      //             << " w: " << int(w) << std::endl;
      // }
      return w;
    }
    }
    return get_max_weight();
  }
  Weight get_min_weight() const { return config_.min_weight; }
  Weight get_max_weight() const { return config_.max_weight; }

private:
  Weight get_weight(double dist_factor) {
    const Weight w = std::min(
        static_cast<Weight>(dist_factor * get_max_weight()), get_max_weight());
    return std::max(get_min_weight(), w);
  }
  WeightGeneratorConfig<Weight> config_;
  std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<Weight> distribution;
};
} // namespace graphs
