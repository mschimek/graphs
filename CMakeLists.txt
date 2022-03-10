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

add_library(base_target INTERFACE)
target_compile_options(base_target INTERFACE -O3 -march=native)
target_link_libraries(base_target INTERFACE MPI::MPI_CXX)
target_include_directories(base_target INTERFACE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include)


target_include_directories(base_target SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/external)

add_library(graphs STATIC
  src/read_from_file.cpp
  src/generate_gnm.cpp
  src/generate_rmat.cpp)
target_link_libraries(graphs PRIVATE base_target kadis generators rmat)
target_include_directories(graphs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(example app/example.cpp)
target_link_libraries(example PRIVATE base_target graphs)

