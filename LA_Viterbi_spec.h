#pragma once

#include "HMM.h"

extern "C" {
    #include <GraphBLAS.h>
}

class LA_Viterbi_spec {
  public:
    explicit LA_Viterbi_spec(const std::string& chmm_file_path);

    std::vector<GrB_Matrix> emit_pr_x_start_pr;
    std::vector<GrB_Matrix> emit_pr_x_trans_pr;
    HMM::Index_t states_num;
    GrB_Matrix result;

    // Writes possible probabilities of seq to result
    void run_Viterbi_spec(const HMM::Emit_vec_t& seq);

    ~LA_Viterbi_spec();
};
