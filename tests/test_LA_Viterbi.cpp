#include "../LA_Viterbi.h"
#include "../HMM_reader.h"

#include <iostream>
int main() {
    auto hmm = read_HMM("../CHMM_data/example.CHMM");
    // Define sequence to examine
    // GGCACTGAA
    // A = 0, C = 1, G = 2, T = 3
    auto seq = HMM::Emit_vec_t{2,2,1,0,1,3,2,0,0};

    auto Viterbi_impl = LA_Viterbi();
    Viterbi_impl.run_Viterbi(hmm, seq);

    // Expected result is
    // (0,0) 25.6574
    // (1,0) 24.4874
    auto prob_0 = HMM::Probability_t(0);
    auto prob_1 = HMM::Probability_t(0);
    GrB_Matrix_extractElement_FP32(&prob_0, Viterbi_impl.result, 0, 0);
    GrB_Matrix_extractElement_FP32(&prob_1, Viterbi_impl.result, 1, 0);
    auto is_test_passed =
        HMM::almost_equal(prob_0, HMM::Probability_t(25.6574)) &&
        HMM::almost_equal(prob_1, HMM::Probability_t(24.4874));
    std::cerr << prob_0 << ' ' << prob_1 << '\n';
    return 1 - static_cast<int>(is_test_passed);
}
