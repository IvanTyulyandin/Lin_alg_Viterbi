#include "../LA_Viterbi_spec.h"
#include "../data_reader.h"

#include <iostream>


int main() {
    auto seq = read_emit_seq("../ess_files/test_seq.ess")[0];

    auto Viterbi_impl = LA_Viterbi_spec("../chmm_files/test_chmm.chmm");
    Viterbi_impl.run_Viterbi_spec(seq);

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
