cmake_minimum_required(VERSION 3.15)
project(Lin_alg_Viterbi)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Wextra -Wpedantic -Werror -Wno-unused-function -Wno-unused-variable)

# Build main part
set(LA_VITERBI_SRC
    HMM.h
    data_reader.h data_reader.cpp
    LA_Viterbi.h LA_Viterbi.cpp)
add_library(LA_Viterbi ${LA_VITERBI_SRC})

find_package(OpenMP REQUIRED)
target_link_libraries(LA_Viterbi graphblas OpenMP::OpenMP_CXX)

# Build executable that can be used for demonstration/benchmarking/etc
add_executable(Lin_alg_Viterbi main.cpp)
target_link_libraries(Lin_alg_Viterbi LA_Viterbi)

# Build tests
enable_testing()
add_subdirectory(tests)

file(GLOB CHMM_DATA ${CMAKE_CURRENT_SOURCE_DIR}/CHMM_data/*.CHMM)
file(COPY ${CHMM_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/CHMM_data/)

file(GLOB ESS_DATA ${CMAKE_CURRENT_SOURCE_DIR}/ess_files/*.ess)
file(COPY ${ESS_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ess_files/)
