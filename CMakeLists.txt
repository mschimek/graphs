cmake_minimum_required(VERSION 3.12)
project(
  graphs
  DESCRIPTION "Unified distributed graph generation and input"
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)

add_subdirectory(external/rmat-fork)
add_subdirectory(external/KaGen)
add_subdirectory(external/KaDiS)

add_library(graphs_base_target INTERFACE)
target_compile_options(graphs_base_target INTERFACE -O3 -march=native -std=c++17)
target_link_libraries(graphs_base_target INTERFACE MPI::MPI_CXX)
target_include_directories(graphs_base_target INTERFACE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include)


target_include_directories(graphs_base_target SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/external)

add_library(graphs STATIC
  src/read_from_file.cpp
  src/generate_gnm.cpp
  src/generate_rgg.cpp
  src/generate_rmat.cpp)
target_link_libraries(graphs PRIVATE graphs_base_target kadis generators rmat)
target_include_directories(graphs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(example app/example.cpp)
add_executable(generate_gnm app/generate_gnm.cpp)
add_executable(generate_rmat app/generate_rmat.cpp)
add_executable(generate_rgg app/generate_rgg.cpp)
target_link_libraries(example PRIVATE graphs_base_target graphs)
target_link_libraries(generate_gnm PRIVATE graphs_base_target graphs)
target_link_libraries(generate_rmat PRIVATE graphs_base_target graphs)
target_link_libraries(generate_rgg PRIVATE graphs_base_target graphs)

