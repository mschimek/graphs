#include "interface.hpp"

#include <cstring>
#include <fstream>
#include <random>
#include <sstream>

namespace graphs {
template <typename T>
void append(std::vector<char>& data, T t, std::size_t num_bytes) {
  const auto prev_size = data.size();
  data.resize(data.size() + num_bytes);
  auto dest = data.data() + prev_size;
  std::memcpy(dest, &t, num_bytes);
}
bool process_chunk(std::ifstream& in, std::ofstream& out,
                   std::size_t chunk_size, Weight max_size,
                   std::size_t num_bytes_vid, std::size_t num_weight_bytes,
                   bool is_graph_directed, bool is_graph_weighted) {
  static std::mt19937 gen(0);
  static std::uniform_int_distribution<std::uint8_t> distrib(1, max_size);
  std::string line;
  std::vector<char> data;
  std::size_t bytes_per_edge = 2 * num_bytes_vid + num_weight_bytes;
  data.reserve(2 * chunk_size * bytes_per_edge);
  bool has_next_line = false;
  for (std::size_t i = 0;
       i < chunk_size &&
       (has_next_line = static_cast<bool>(std::getline(in, line, '\n')));
       ++i) {
    std::stringstream ss(line);
    VId src, dst;
    std::uint64_t weight = distrib(gen);
    ss >> src >> dst;
    if (is_graph_weighted)
      ss >> weight;
    std::cout << src << " " << dst << " " << weight << std::endl;

    if (src == dst) {
      continue;
    }
    if (!is_graph_directed && !is_graph_weighted && src > dst) {
      continue;
    }
    append(data, src, num_bytes_vid);
    append(data, dst, num_bytes_vid);
    append(data, weight, num_weight_bytes);
    if (!is_graph_weighted) {
      append(data, dst, num_bytes_vid);
      append(data, src, num_bytes_vid);
      append(data, weight, num_weight_bytes);
    }
  }
  out.write(data.data(), data.size());
  return has_next_line;
}


std::size_t get_file_size(std::string filename) {
  std::ifstream end_of_file(filename.c_str(), std::ios::binary | std::ios::ate);
  std::fstream::pos_type file_size;
  if (end_of_file) {
    return end_of_file.tellg();
  }
  return 0;
}

std::vector<WEdge> get_edges(const Configs& config) {
  const std::size_t num_input_bytes = get_file_size(config.infile);
  std::ifstream in{config.infile.c_str(), std::ios::in};
  std::vector<unsigned char> data(num_input_bytes);
  std::vector<WEdge> edges;
  const std::size_t edge_size =
      2 * config.num_vid_bytes_in + config.num_weight_bytes_out;
  edges.reserve(num_input_bytes / edge_size);
  for (std::size_t i = 0; i < data.size(); i += edge_size) {
    auto it = data.data() + i;
    VId src = 0;
    VId dst = 0;
    uint8_t w = 0;
    std::memcpy(&src, it, config.num_vid_bytes_in);
    it += config.num_vid_bytes_in;
    std::memcpy(&dst, it, config.num_vid_bytes_in);
    it += config.num_vid_bytes_in;
    std::memcpy(&w, it, config.num_weight_bytes_in);
    edges.emplace_back(src, dst, w);
  }
  return edges;
}
void write_edges(const Configs& config, const std::vector<WEdge>& edges) {
  std::size_t edge_size_byte =
      2 * config.num_vid_bytes_out + config.num_weight_bytes_out;
  std::vector<char> data;
  data.reserve(edge_size_byte);
  auto it = data.begin();
  for (const auto& edge : edges) {
    append(data, edge.get_src(), config.num_vid_bytes_out);
    append(data, edge.get_dst(), config.num_vid_bytes_out);
    append(data, edge.get_weight(), config.num_weight_bytes_out);
  }
  std::ofstream out = std::ofstream(
      config.outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  out.write(data.data(), data.size());
}

void read_sort_write(const Configs& config) {
  std::ifstream in{config.infile.c_str(), std::ios::in};
  std::ofstream out = std::ofstream(
      config.outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  auto edges = get_edges(config);
  ips4o::parallel::sort(edges.begin(), edges.end(), SrcDstWeightOrder<WEdge>{});
  write_edges(config, edges);
}

void add_weight_and_back_edges_in_directed_graph(const std::string& infile,
                                                 const std::string& outfile,
                                                 Weight max_weight,
                                                 std::size_t num_VId_bytes,
                                                 std::size_t num_weight_bytes) {
  std::ifstream in{infile.c_str(), std::ios::in};
  std::ofstream out = std::ofstream(
      outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  const bool is_graph_directed = true;
  const bool is_graph_weighted = false;
  while (process_chunk(in, out, 10'000'000, max_weight, num_VId_bytes,
                       num_weight_bytes, is_graph_directed, is_graph_weighted))
    ;
}
void add_weight_and_back_edges_in_undirected_graph(
    const std::string& infile, const std::string& outfile, Weight max_weight,
    std::size_t num_VId_bytes, std::size_t num_weight_bytes) {
  std::ifstream in{infile.c_str(), std::ios::in};
  std::ofstream out = std::ofstream(
      outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  const bool is_graph_directed = false;
  const bool is_graph_weighted = false;
  while (process_chunk(in, out, 10'000'000, max_weight, num_VId_bytes,
                       num_weight_bytes, is_graph_directed, is_graph_weighted))
    ;
}

void transform_undirected_graph_to_sorted_binary(const std::string& infile,
                                                 const std::string& outfile,
                                                 Weight max_weight,
                                                 std::size_t num_VId_bytes,
                                                 std::size_t num_weight_bytes) {
  std::ifstream in{infile.c_str(), std::ios::in};
  std::ofstream out = std::ofstream(
      outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  const bool is_graph_directed = false;
  const bool is_graph_weighted = true;
  while (process_chunk(in, out, 10'000'000, max_weight, num_VId_bytes,
                       num_weight_bytes, is_graph_directed, is_graph_weighted))
    ;
}

inline VId get_number(char* ptr) {
  static_assert(sizeof(VId) == sizeof(uint64_t));
  VId read = -1;
  std::memcpy(&read, ptr, sizeof(VId));
  return le64toh(read);
}

template <typename Container>
Container read_weighted_binary(const std::string& infile,
                               std::size_t num_VId_bytes,
                               std::size_t num_weight_bytes) {
  using EdgeType = typename Container::value_type;
  std::ifstream in{infile.c_str(), std::ios::binary | std::ios::in};
  std::ifstream end_of_file(infile.c_str(), std::ios::binary | std::ios::ate);
  std::fstream::pos_type file_size;
  if (end_of_file) {
    file_size = end_of_file.tellg();
  }
  std::vector<char> data(file_size);
  in.read(reinterpret_cast<char*>(data.data()), file_size);
  std::cout << "file size: " << file_size << " "
            << sizeof(std::fstream::pos_type) << std::endl;
  std::size_t step_size =
      2 * num_VId_bytes + num_weight_bytes; // src, dst, weight
  for (std::size_t i = 0; i < data.size(); i += step_size) {
    auto it = data.data() + i;
    VId src = 0;
    VId dst = 0;
    uint8_t w = 0;
    std::memcpy(&src, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&dst, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&w, it, num_weight_bytes);
    WEdge e{src, dst, w};
    std::cout << e << std::endl;
  }
}
} // namespace graphs
