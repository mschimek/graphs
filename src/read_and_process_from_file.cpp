#include "interface.hpp"

#include <cstring>
#include <fstream>
#include <random>
#include <sstream>

namespace graphs {
template <typename T> void append(std::vector<char>& data, T t, std::size_t num_bytes) {
  const auto prev_size = data.size();
  data.resize(data.size() + num_bytes);
  auto dest = data.data() + prev_size;
  std::memcpy(dest, &t, num_bytes);
}
bool process_chunk(std::ifstream& in, std::ofstream& out,
                   std::size_t chunk_size, Weight max_size, std::size_t num_bytes_vid) {
  static std::mt19937 gen(0);
  static std::uniform_int_distribution<std::uint8_t> distrib(1, max_size);
  std::string line;
  std::vector<char> data;
  data.reserve(chunk_size * 2 * 13);
  bool has_next_line = false;
  for (std::size_t i = 0;
       i < chunk_size &&
       (has_next_line = static_cast<bool>(std::getline(in, line, '\n')));
       ++i) {
    std::stringstream ss(line);
    VId src, dst;
    ss >> src >> dst;
    if (src == dst) {
      continue;
    }
    std::uint8_t w = distrib(gen);
    append(data, src, num_bytes_vid);
    append(data, dst, num_bytes_vid);
    append(data, w, 1);
    append(data, dst, num_bytes_vid);
    append(data, src, num_bytes_vid);
    append(data, w, 1);
  }
  out.write(data.data(), data.size());
  return has_next_line;
}
void add_weight_and_back_edges(const std::string& infile,
                               const std::string& outfile, Weight max_weight,
                               std::size_t num_VId_bytes) {
  std::ifstream in{infile.c_str(), std::ios::in};
  std::ofstream out = std::ofstream(
      outfile.c_str(), std::ios::out | std::ios::binary | std::ios::app);
  while (process_chunk(in, out, 10'000'000, max_weight, num_VId_bytes))
    ;
}

inline VId get_number(char* ptr) {
  static_assert(sizeof(VId) == sizeof(uint64_t));
  VId read = -1;
  std::memcpy(&read, ptr, sizeof(VId));
  return le64toh(read);
}

template<typename Container>
Container read_weighted_binary(const std::string& infile,
                          std::size_t num_VId_bytes) {
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
  std::size_t step_size = 2 * num_VId_bytes + 1; // src, dst, weight
  for (std::size_t i = 0; i < data.size(); i += step_size) {
    auto it = data.data() + i;
    VId src = 0;
    VId dst = 0;
    uint8_t w = 0;
    std::memcpy(&src, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&dst, it, num_VId_bytes);
    it += num_VId_bytes;
    std::memcpy(&w, it, 1);
    WEdge e{src, dst, w};
    std::cout << e << std::endl;
  }
}
} // namespace graphs
