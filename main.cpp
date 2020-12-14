#include "benchmark_helper.h"
#include <iostream>


int main() {

    launch_GraphBLAS();

    auto ess = read_emit_seq("./ess_files/emit_3_3500_20.ess");
    auto file_with_chmms = std::string("./chmm_files/");
    constexpr int TIMES_TO_RUN = 1;

    // Benchmark non-specialized and specialized versions
    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algorithm_selector::LA>(
        file_with_chmms, ess);

    benchmark_with_chmms_in_folder<TIMES_TO_RUN, Algorithm_selector::LA_spec>(
        file_with_chmms, ess);

    stop_GraphBLAS();

    return 0;
}
