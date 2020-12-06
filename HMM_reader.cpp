#include "HMM_reader.h"

#include <cmath>
#include <iostream>
#include <fstream>


// HMM format
//
// HMM states amount N
// N starting probabilities at one string
// emit symbols amount M
// N strings with M probabilities
// transitions amount T
// T lines of: src_state dest_state probability


HMM read_HMM(const std::string& HMM_file_name) {
    auto file = std::ifstream(HMM_file_name);
    if (file.fail()) {
        std::cerr << "Failed to open file with HMM: " << HMM_file_name << '\n';
        return HMM{};
    }

    auto hmm = HMM{};
    auto prob_from_file = HMM::Probability_t(0);

    // Read number of states and probabilities to be start state
    file >> hmm.states_num;
    hmm.start_probabilities.reserve(hmm.states_num);

    for (size_t i = 0; i < hmm.states_num; ++i) {
        file >> prob_from_file;
        prob_from_file = HMM::to_neg_log(prob_from_file);
        hmm.start_probabilities.push_back(prob_from_file);
    }

    // Read info about emission symbols
    file >> hmm.emit_num;
    hmm.emissions.reserve(hmm.emit_num * hmm.states_num);

    for (size_t i = 0; i < hmm.emit_num * hmm.states_num; ++i) {
        file >> prob_from_file;
        prob_from_file = HMM::to_neg_log(prob_from_file);
        hmm.emissions.push_back(prob_from_file);
    }

    // Read graph edges info as triples:
    // source_state destination_state transition_probability
    file >> hmm.trans_num;
    hmm.trans_rows = HMM::Index_vec_t();
    hmm.trans_rows.reserve(hmm.trans_num);
    hmm.trans_cols = HMM::Index_vec_t();
    hmm.trans_cols.reserve(hmm.trans_num);

    auto src = size_t(0);
    auto dst = size_t(0);
    for (size_t i = 0; i < hmm.trans_num; ++i) {
        file >> src >> dst >> prob_from_file;
        prob_from_file = HMM::to_neg_log(prob_from_file);
        hmm.trans_rows.push_back(src);
        hmm.trans_cols.push_back(dst);
        hmm.trans_probs.push_back(prob_from_file);
    }

    file.close();
    return hmm;
}
