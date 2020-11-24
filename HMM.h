#include <vector>

class HMM {
  public:
    using Probability_t = float;

    // encode states and emit symbols as numbers
    std::vector<Probability_t> transitions;
    std::vector<Probability_t> emissions;
    std::vector<Probability_t> start_probabilities;
};
