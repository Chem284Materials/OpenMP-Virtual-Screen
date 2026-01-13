cmake_minimum_required(VERSION 3.14)

# Project name and version
project(OpenMPVirtualScreen VERSION 1.0
	    DESCRIPTION "Parallelize a VS with OpenMP"
		LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Add the executables
add_executable(vector_norm src/vector_norm.cpp)
add_executable(serial_docking src/serial_docking.cpp)
add_executable(parallel_docking_manual_reduction src/parallel_docking_manual_reduction.cpp)
add_executable(parallel_docking_custom_reduction src/parallel_docking_custom_reduction.cpp)
add_executable(simd_parallel_docking src/simd_parallel_docking.cpp)

# Add include for serial_docking
target_include_directories(serial_docking PRIVATE include/)
target_include_directories(parallel_docking_manual_reduction PRIVATE include/)
target_include_directories(parallel_docking_custom_reduction PRIVATE include/)
target_include_directories(simd_parallel_docking PRIVATE include/)

# Find OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(vector_norm PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(parallel_docking_manual_reduction PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(parallel_docking_custom_reduction PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(simd_parallel_docking PUBLIC OpenMP::OpenMP_CXX)

target_compile_options(simd_parallel_docking PRIVATE -O3 -mavx2 -fopt-info-vec-optimized)
