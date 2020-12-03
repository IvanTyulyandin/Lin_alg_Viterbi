#pragma once

#include <cmath>
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

    // Functions to work with Probability_t

    static bool almost_equal(HMM::Probability_t x, HMM::Probability_t y) {
        return std::fabs(x - y) <= 0.0001;
    }

    static HMM::Probability_t to_neg_log(HMM::Probability_t x) {
        return -1 * std::log2(x);
    }
};
