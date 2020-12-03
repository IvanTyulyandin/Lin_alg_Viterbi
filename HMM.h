#pragma once

#include <vector>

class HMM {
  public:
    using Probability_t = float;
    using Index_t = size_t;
    using Emit_t = size_t;

    Index_t states_num;
    Index_t emit_num;
    Index_t trans_num;
    // encode states and emit symbols as numbers
    std::vector<Index_t> trans_rows;
    std::vector<Index_t> trans_cols;
    std::vector<Probability_t> trans_probs;
    std::vector<Probability_t> emissions;
    std::vector<Probability_t> start_probabilities;
};
