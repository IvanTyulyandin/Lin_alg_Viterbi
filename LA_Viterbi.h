#pragma once

#include "HMM.h"

extern "C" {
    #include <GraphBLAS.h>
}

class LA_Viterbi {
  public:
    LA_Viterbi();

    GrB_Matrix result;

    // Writes possible probabilities of seq to result
    void run_Viterbi(HMM& hmm, HMM::Emit_vec_t seq);

    ~LA_Viterbi();
};
