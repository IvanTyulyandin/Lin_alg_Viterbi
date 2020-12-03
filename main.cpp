#include "HMM_reader.h"
#include "LA_Viterbi.h"

#include <iostream>


int main() {
    auto hmm = read_HMM("example.CHMM");
    // Define sequence to examine
    // GGCACTGAA
    // A = 0, C = 1, G = 2, T = 3
    auto seq = std::vector<HMM::Emit_t>{2,2,1,0,1,3,2,0,0};
    // Expected result is
    // (0,0) 25.6574
    // (1,0) 24.4874
    auto Viterbi_impl = LA_Viterbi();
    Viterbi_impl.run_Viterbi(hmm, seq);
    auto prob = HMM::Probability_t(0);
    GrB_Matrix_extractElement_FP32(&prob, Viterbi_impl.result,  0, 0);
    if (!HMM::almost_equal(prob, float(25.6574))) {
        std::cerr << "Wrong (0,0) or result probabilities!\n"
            << "(0,0) is " << prob << '\n';
    }
    GrB_Matrix_extractElement_FP32(&prob, Viterbi_impl.result, 1, 0);
    if (!HMM::almost_equal(prob, float(24.4874))) {
        std::cerr << "Wrong (0,1) or result probabilities!\n"
            << "(1,0) is " << prob << '\n';
    }
    return 0;
}
