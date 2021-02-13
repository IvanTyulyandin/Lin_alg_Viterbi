#include "benchmark_helper.h"
#include <iostream>

constexpr int TIMES_TO_RUN = 100;

void run_benchmark(const std::string& file_with_chmms, const HMM::Seq_vec_t& ess) {

    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algo_selector::LA>(file_with_chmms, ess);

    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algo_selector::LA_spec>(file_with_chmms, ess);

    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algo_selector::LA_assoc_spec, 2>(file_with_chmms,
                                                                                  ess);

    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algo_selector::LA_assoc_spec, 3>(file_with_chmms,
                                                                                  ess);
}

int main() {

    launch_GraphBLAS();

    auto ess_files = std::vector{"./ess_files/emit_3_3500_20.ess", "./ess_files/emit_3_7000_20.ess",
                                 "./ess_files/covid-19.ess"};

    for (auto ess_file : ess_files) {
        std::cout << "----------------------------\n"
                  << ess_file << "\n----------------------------\n";
        auto ess = read_emit_seq(ess_file);
        auto file_with_chmms = std::string("./chmm_files/");
        run_benchmark(file_with_chmms, ess);
    }

    stop_GraphBLAS();

    return 0;
}
