#pragma once

#include "GraphBLAS_manager.h"
#include "HMM.h"

class LA_Viterbi {
  public:
    LA_Viterbi();

    GrB_Matrix result;

    // Writes possible probabilities of seq to result
    void run_Viterbi(HMM& hmm, const HMM::Emit_vec_t& seq);

    ~LA_Viterbi();
};
