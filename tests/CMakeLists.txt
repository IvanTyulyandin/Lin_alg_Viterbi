cmake_minimum_required(VERSION 3.15)
project(tests)

add_executable(test_chmm_reader test_chmm_reader.cpp)
target_link_libraries(test_chmm_reader LA_Viterbi)

add_test(test_chmm_reader test_chmm_reader)

add_executable(test_ess_reader test_ess_reader.cpp)
target_link_libraries(test_ess_reader LA_Viterbi)

add_test(test_ess_reader test_ess_reader)

add_executable(test_LA_Viterbi test_LA_Viterbi.cpp)
target_link_libraries(test_LA_Viterbi LA_Viterbi)

add_test(test_LA_Viterbi test_LA_Viterbi)

add_executable(test_LA_Viterbi_spec test_LA_Viterbi_spec.cpp)
target_link_libraries(test_LA_Viterbi_spec LA_Viterbi)

add_test(test_LA_Viterbi_spec test_LA_Viterbi_spec)
