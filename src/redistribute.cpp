#include "include/utils.hpp"
#include "interface.hpp"

namespace graphs {

std::vector<WEdge> alltoall(const std::vector<WEdge>& edges,
                            const std::vector<std::size_t>& send_counts_,
                            MPIComm comm) {
  const int size = comm.size;

  std::vector<WEdge> recv_buf;
  std::vector<int> send_counts(size);
  std::vector<int> recv_counts(size);
  std::vector<int> send_displs(size);
  std::vector<int> recv_displs(size);
  for (size_t i = 0; i < send_counts.size(); ++i) {
    send_counts[i] = send_counts_[i];
  }

  std::exclusive_scan(send_counts.begin(), send_counts.end(),
                      send_displs.begin(), 0);
  const std::size_t total_send_count = send_displs.back() + send_counts.back();
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               comm.comm);
  std::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                      recv_displs.begin(), 0);
  const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();

  recv_buf.resize(total_recv_count);
  const auto wedge_mpi_type = WEdge::MPI_Type{};
  MPI_Alltoallv(edges.data(), send_counts.data(), send_displs.data(),
                wedge_mpi_type.get_mpi_type(), recv_buf.data(),
                recv_counts.data(), recv_displs.data(),
                wedge_mpi_type.get_mpi_type(), comm.comm);
  return recv_buf;
}
std::vector<WEdge> redistribute(const std::vector<WEdge>& local_data,
                                std::uint64_t num_edges, MPIComm comm) {
  std::vector<std::uint64_t> num_requested_elements(comm.size, 0);
  std::vector<std::uint64_t> num_actual_elements(comm.size, 0);
  std::vector<std::uint64_t> in_scan_requested_elements(comm.size, 0);
  std::vector<std::uint64_t> ex_scan_actual_elements(comm.size, 0);
  const std::uint64_t actual_num_edges = local_data.size();
  MPI_Allgather(&num_edges, 1, MPI_UINT64_T, num_requested_elements.data(), 1,
                MPI_UINT64_T, comm.comm);
  MPI_Allgather(&actual_num_edges, 1, MPI_UINT64_T, num_actual_elements.data(),
                1, MPI_UINT64_T, comm.comm);
  std::inclusive_scan(num_requested_elements.begin(),
                      num_requested_elements.end(),
                      in_scan_requested_elements.begin());
  std::exclusive_scan(num_actual_elements.begin(), num_actual_elements.end(),
                      ex_scan_actual_elements.begin(), 0ull);

  std::size_t num_remaining_elements = local_data.size();
  std::size_t offset = ex_scan_actual_elements[comm.rank];
  std::vector<std::uint64_t> send_counts(comm.size, 0);
  for (std::size_t i = 0; i < comm.size && num_remaining_elements > 0; ++i) {
    if (in_scan_requested_elements[i] <=
        offset) { // could also use a binary search;
      continue;
    }
    const std::size_t max_elements_cur_rank =
        in_scan_requested_elements[i] - offset;
    const std::size_t elements_to_cur_rank =
        std::min(num_remaining_elements, max_elements_cur_rank);
    num_remaining_elements -= elements_to_cur_rank;
    offset += elements_to_cur_rank;
    send_counts[i] = elements_to_cur_rank;
  }
  //graphs::execute_in_order(comm, [&]() {
  //  std::cout << "num_edges: " << num_edges
  //            << " actual_num_edges: " << actual_num_edges << std::endl;
  //  for (std::size_t i = 0; i < comm.size; ++i) {
  //    std::cout << "req: " << num_requested_elements[i]
  //              << " actual: " << num_actual_elements[i]
  //              << " in scan requested: " << in_scan_requested_elements[i]
  //              << " ex scan actual: " << ex_scan_actual_elements[i] << " "
  //              << std::endl;
  //  }
  //});
  auto res = alltoall(local_data, send_counts, comm);
  if(res.size() != num_edges) {
    std::cout << "wrong redistribution" << std::endl;
    std::abort();
  }
  return res;
}
} // namespace graphs
