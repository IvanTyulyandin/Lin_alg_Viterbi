#include "../data_reader.h"

int main() {
    auto emitted_sequences = read_emit_seq("../ess_files/test_seq.ess");
    auto is_test_passed = emitted_sequences.size() == 2 &&
                          emitted_sequences[0] == HMM::Emit_vec_t{2, 2, 1, 0, 1, 3, 2, 0, 0} &&
                          emitted_sequences[1] == HMM::Emit_vec_t{3, 2, 1, 0};

    return 1 - static_cast<int>(is_test_passed);
}
