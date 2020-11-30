#include "HMM_reader.h"

#include <iostream>
#include <limits>
#include <cmath>
#include <utility> // std::swap
extern "C" {
    #include <GraphBLAS.h>
}

void check_for_error(const GrB_Info& info) {
    if (! (info == GrB_Info::GrB_SUCCESS || info == GrB_Info::GrB_NO_VALUE)) {
        printf("info: %d error: %s\n", info, GrB_error());
    }
}


int main() {
    auto hmm = read_HMM("example.CHMM");

    // Init GraphBLAS
    auto info = GrB_init(GrB_Mode::GrB_NONBLOCKING);
    check_for_error(info);

    // Define GraphBLAS matrices

    // Start/current probabilities
    auto cur_probs = GrB_Matrix();
    info = GrB_Matrix_new(&cur_probs, GrB_FP32, hmm.states_num, 1);
    check_for_error(info);

    auto cur_probs_rows_ind = std::vector<GrB_Index>(hmm.states_num);
    auto cur_probs_cols_ind = std::vector<GrB_Index>(hmm.states_num);
    for (size_t i = 0; i < hmm.states_num; ++i) {
        cur_probs_rows_ind[i] = i;
        cur_probs_cols_ind[i] = 0;
    }
    info = GrB_Matrix_build_FP32(cur_probs,
        cur_probs_rows_ind.data(), cur_probs_cols_ind.data(), hmm.start_probabilities.data(),
        hmm.states_num,
        GrB_FIRST_FP32);
    check_for_error(info);

    // Emission probabilities matrices
    auto em_probs = std::vector<GrB_Matrix>(hmm.emit_num);
    auto emit_ind = std::vector<GrB_Index>(hmm.states_num);
    auto emit_data = std::vector<HMM::Probability_t>(hmm.states_num);

    // Diagonal matrix indicies
    for (size_t i = 0; i < hmm.states_num; ++i) {
        emit_ind[i] = i;
    }

    for (size_t i = 0; i < hmm.emit_num; ++i) {
        auto& m = em_probs[i];
        m = GrB_Matrix();
        info = GrB_Matrix_new(&m, GrB_FP32, hmm.states_num, hmm.states_num);
        check_for_error(info);

        auto offset = i;
        for (size_t j = 0; j < hmm.states_num; ++j) {
            emit_data[j] = hmm.emissions[offset];
            offset += hmm.emit_num;
        }

        info = GrB_Matrix_build_FP32(
            m,
            emit_ind.data(), emit_ind.data(), emit_data.data(),
            hmm.states_num,
            GrB_FIRST_FP32);
        check_for_error(info);
    }

    // Transposed transition matrix
    std::swap(hmm.trans_rows, hmm.trans_cols);

    auto transitions = GrB_Matrix();
    info = GrB_Matrix_new(&transitions, GrB_FP32, 2, 2);
    check_for_error(info);

    info = GrB_Matrix_build_FP32(
        transitions,
        hmm.trans_rows.data(),hmm.trans_cols.data(), hmm.trans_probs.data(),
        hmm.trans_num,
        GrB_FIRST_FP32);
    check_for_error(info);

    // Define sequence to examine
    // GGCACTGAA
    // A = 0, C = 1, G = 2, T = 3
    auto seq = std::vector<size_t>{2,2,1,0,1,3,2,0,0};

    auto prob_x_trans = GrB_Matrix();
    info = GrB_Matrix_new(&prob_x_trans, GrB_FP32, hmm.states_num, hmm.states_num);
    check_for_error(info);

    auto next_probabilites = GrB_Matrix();
    info = GrB_Matrix_new(&next_probabilites, GrB_FP32, hmm.states_num, 1);
    check_for_error(info);

    // Viterbi algorithm

    // Count emissions for first symbol and start probabilities
    info = GrB_mxm(
        next_probabilites, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
        em_probs[seq[0]], cur_probs, GrB_NULL);
    check_for_error(info);
    std::swap(cur_probs, next_probabilites);

    for (size_t i = 1; i < seq.size(); ++i) {
        info = GrB_mxm(
            prob_x_trans, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
            em_probs[seq[i]], transitions, GrB_NULL);
        check_for_error(info);

        info = GrB_mxm(
            next_probabilites, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
            prob_x_trans, cur_probs, GrB_NULL);
        check_for_error(info);

        std::swap(cur_probs, next_probabilites);
    }

    // Print matrix
    // SuiteSPARSE-specific (since GxB, not GrB)
    // Expected result is
    // (0,0) 25.65
    // (1,0) 24.49
    info = GxB_Matrix_fprint(cur_probs, "current_probabilities", GxB_COMPLETE, stdout);
    check_for_error(info);

    // Free resources
    for (auto& m : em_probs) {
        GrB_Matrix_free(&m);
    }
    GrB_Matrix_free(&transitions);
    GrB_Matrix_free(&prob_x_trans);
    GrB_Matrix_free(&cur_probs);
    GrB_Matrix_free(&next_probabilites);

    // Finalize GraphBLAS
    info = GrB_finalize();
    check_for_error(info);

    return 0;
}
