cmake_minimum_required(VERSION 3.15)
project(Lin_alg_Viterbi)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Wextra -Wpedantic -Werror -Wno-unused-function -Wno-unused-variable)

find_package(OpenMP REQUIRED)

add_executable(Lin_alg_Viterbi main.cpp HMM.h HMM_reader.h HMM_reader.cpp)
target_link_libraries(Lin_alg_Viterbi graphblas OpenMP::OpenMP_CXX)

file(GLOB TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/*.HMM)
file(COPY ${TEST_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
