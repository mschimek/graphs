cmake_minimum_required(VERSION 3.12)
project(
  graphs
  DESCRIPTION "Unified distributed graph generation and input"
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

add_subdirectory(external/rmat-fork)
add_subdirectory(external/KaDiS)
add_subdirectory(external/KAGEN)

add_library(graphs_base_target INTERFACE)
target_compile_options(graphs_base_target INTERFACE -O3 -march=native -std=c++17)
target_link_libraries(graphs_base_target INTERFACE MPI::MPI_CXX OpenMP::OpenMP_CXX KaGen::KaGen)
target_include_directories(graphs_base_target INTERFACE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include)


target_include_directories(graphs_base_target SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/external)

add_library(graphs STATIC
  src/read_from_file.cpp
  src/generate_gnm.cpp
  src/generate_rgg.cpp
  src/generate_rhg.cpp
  src/generate_rmat.cpp
  src/generate_grid.cpp
  src/redistribute.cpp
  src/read_and_process_from_file.cpp
  src/utils.cpp)
target_link_libraries(graphs PRIVATE graphs_base_target  rmat kadis)
target_include_directories(graphs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(example app/example.cpp)
add_executable(generate_gnm app/generate_gnm.cpp)
add_executable(generate_rmat app/generate_rmat.cpp)
add_executable(generate_rgg app/generate_rgg.cpp)
add_executable(generate_rhg app/generate_rhg.cpp)
add_executable(generate_grid app/generate_grid.cpp)
add_executable(generate_grid3d app/generate_grid3d.cpp)
add_executable(process_file app/process_file.cpp)
add_executable(read_weighted_binary_file app/read_weighted_binary_file.cpp)
add_executable(read_sort_write app/read_sort_write.cpp)
add_executable(redistribute app/redistribute.cpp)
target_link_libraries(example PRIVATE graphs_base_target graphs)
target_link_libraries(generate_gnm    PRIVATE graphs_base_target graphs)
target_link_libraries(generate_rmat   PRIVATE graphs_base_target graphs)
target_link_libraries(generate_rgg    PRIVATE graphs_base_target graphs)
target_link_libraries(generate_rhg    PRIVATE graphs_base_target graphs)
target_link_libraries(generate_grid   PRIVATE graphs_base_target graphs)
target_link_libraries(generate_grid3d PRIVATE graphs_base_target graphs)
target_link_libraries(redistribute PRIVATE graphs_base_target graphs)
target_link_libraries(process_file PRIVATE graphs_base_target graphs tlx kadis) # depends on tlx in KaDiS ... should be fixed
target_link_libraries(read_weighted_binary_file PRIVATE graphs_base_target graphs tlx kadis) # depends on tlx in KaDiS ... should be fixed
target_link_libraries(read_sort_write PRIVATE graphs_base_target graphs tlx kadis) # depends on tlx in KaDiS ... should be fixed

