cmake_minimum_required(VERSION 3.15)
project(tests)

add_executable(test_HMM_reader test_HMM_reader.cpp)
target_link_libraries(test_HMM_reader LA_Viterbi)

add_test(test_HMM_reader test_HMM_reader)
