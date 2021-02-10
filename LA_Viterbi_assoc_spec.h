#pragma once

#include "GraphBLAS_manager.h"
#include "HMM.h"

#include <unordered_map>

class LA_Viterbi_assoc_spec {
  public:
    // "level" is the maximum depth for precalc_observation_handlers
    // if level == 1, it makes no sence, same matrices are in emit_pr_x_trans_pr
    // if level == 2, it will save all possible (emit_x_tr * emit_x_tr)
    // and so on
    explicit LA_Viterbi_assoc_spec(const std::string& chmm_file_path, size_t level);

    using Obs_handler_t = std::unordered_map<std::vector<size_t>, GrB_Matrix, HMM::Emit_vec_hasher>;

    std::vector<GrB_Matrix> emit_pr_x_start_pr;
    std::vector<GrB_Matrix> emit_pr_x_trans_pr;
    Obs_handler_t precalc_obs_handlers;
    HMM::Index_t states_num;
    size_t level;
    GrB_Matrix result;

    // Writes possible probabilities of seq to result
    void run_Viterbi_assoc_spec(const HMM::Emit_vec_t& seq);

    ~LA_Viterbi_assoc_spec();
};
