cmake_minimum_required(VERSION 3.15)
project(tests)

add_executable(test_HMM_reader test_HMM_reader.cpp)
target_link_libraries(test_HMM_reader LA_Viterbi)

add_test(test_HMM_reader test_HMM_reader)

add_executable(test_LA_Viterbi test_LA_Viterbi.cpp)
target_link_libraries(test_LA_Viterbi LA_Viterbi)

add_test(test_LA_Viterbi test_LA_Viterbi)

add_executable(test_ess_reader test_ess_reader.cpp)
target_link_libraries(test_ess_reader LA_Viterbi)

add_test(test_ess_reader test_ess_reader)
