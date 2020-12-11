#include "LA_Viterbi_spec.h"

#include "data_reader.h"

namespace {
    void check_for_error(const GrB_Info& info) {
        if (! (info == GrB_Info::GrB_SUCCESS || info == GrB_Info::GrB_NO_VALUE)) {
            printf("info: %d error: %s\n", info, GrB_error());
        }
    }
}


LA_Viterbi_spec::LA_Viterbi_spec(const std::string& chmm_file_path) {
    // Init GraphBLAS
    auto info = GrB_init(GrB_Mode::GrB_NONBLOCKING);
    check_for_error(info);

    auto hmm = read_HMM(chmm_file_path);
    states_num = hmm.states_num;
    result = GrB_Matrix();
    info = GrB_Matrix_new(&result, GrB_FP32, states_num, 1);
    check_for_error(info);

    // Transposed transition matrix
    (hmm.trans_rows).swap(hmm.trans_cols);

    auto transitions = GrB_Matrix();
    info = GrB_Matrix_new(&transitions, GrB_FP32, 2, 2);
    check_for_error(info);

    info = GrB_Matrix_build_FP32(
        transitions,
        hmm.trans_rows.data(),hmm.trans_cols.data(), hmm.trans_probs.data(),
        hmm.trans_num,
        GrB_FIRST_FP32);
    check_for_error(info);

    // Read info about states with
    // non zero probabilities to be start

    auto n_zeroes_ind = std::vector<GrB_Index>(states_num, 0);
    auto from_0_to_n_ind = std::vector<GrB_Index>(states_num);
    for (size_t i = 0; i < states_num; ++i) {
        from_0_to_n_ind[i] = i;
    }

    auto start_probs = GrB_Matrix();
    info = GrB_Matrix_new(&start_probs, GrB_FP32, states_num, 1);
    info = GrB_Matrix_build_FP32(start_probs,
        from_0_to_n_ind.data(), n_zeroes_ind.data(), hmm.start_probabilities.data(),
        hmm.states_num,
        GrB_FIRST_FP32);
    check_for_error(info);

    emit_pr_x_start_pr = std::vector<GrB_Matrix>(hmm.emit_num);
    emit_pr_x_trans_pr = std::vector<GrB_Matrix>(hmm.emit_num);
    for (size_t i = 0; i < hmm.emit_num; ++i) {
        info = GrB_Matrix_new(
            &(emit_pr_x_start_pr[i]), GrB_FP32, states_num, 1);
        check_for_error(info);
        info = GrB_Matrix_new(
            &(emit_pr_x_trans_pr[i]), GrB_FP32, states_num, states_num);
        check_for_error(info);
    }

    auto emit_data = HMM::Prob_vec_t(states_num);
    auto emit_probs_diag_mat = GrB_Matrix();
    info = GrB_Matrix_new(&emit_probs_diag_mat, GrB_FP32, states_num, states_num);

    for (size_t i = 0; i < hmm.emit_num; ++i) {

        auto offset = i;
        for (size_t j = 0; j < states_num; ++j) {
            emit_data[j] = hmm.emissions[offset];
            offset += hmm.emit_num;
        }

        info = GrB_Matrix_build_FP32(
            emit_probs_diag_mat,
            from_0_to_n_ind.data(), from_0_to_n_ind.data(), emit_data.data(),
            hmm.states_num,
            GrB_FIRST_FP32);
        check_for_error(info);

        info = GrB_mxm(
            emit_pr_x_start_pr[i], GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
            emit_probs_diag_mat, start_probs, GrB_NULL);
        check_for_error(info);

        info = GrB_mxm(
            emit_pr_x_trans_pr[i], GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
            emit_probs_diag_mat, transitions, GrB_NULL);
        check_for_error(info);

        info = GrB_Matrix_clear(emit_probs_diag_mat);
        check_for_error(info);
    }

    GrB_Matrix_free(&transitions);
    GrB_Matrix_free(&emit_probs_diag_mat);
}


void LA_Viterbi_spec::run_Viterbi_spec(const HMM::Emit_vec_t& seq) {
    auto info = GrB_Matrix_clear(result);
    check_for_error(info);

    // Start Viterbi algorithm for seq[0]
    info = GrB_Matrix_dup(&result, emit_pr_x_start_pr[seq[0]]);
    check_for_error(info);

    auto next_probs = GrB_Matrix();
    info = GrB_Matrix_new(&next_probs, GrB_FP32, states_num, 1);

    // Viterbi algorithm for the rest of seq
    for (size_t i = 1; i < seq.size(); ++i) {
        info = GrB_mxm(
            next_probs, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP32,
            emit_pr_x_trans_pr[seq[i]], result, GrB_NULL);
        check_for_error(info);

        std::swap(next_probs, result);
    }

    GrB_Matrix_free(&next_probs);

    return;
}


LA_Viterbi_spec::~LA_Viterbi_spec() {
    // Free internal matrices
    for (auto& m: emit_pr_x_trans_pr) {
        GrB_Matrix_free(&m);
    }
    for (auto& m: emit_pr_x_start_pr) {
        GrB_Matrix_free(&m);
    }
    GrB_Matrix_free(&result);

    // Finalize GraphBLAS
    auto info = GrB_finalize();
    check_for_error(info);
}
