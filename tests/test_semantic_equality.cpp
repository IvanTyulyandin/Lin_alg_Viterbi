#include "../LA_Viterbi.h"
#include "../LA_Viterbi_spec.h"
#include "../data_reader.h"

#include <experimental/filesystem>
#include <iostream>

namespace {
constexpr auto chmm_folder = "../chmm_files/";
constexpr auto seq_file_name = "../ess_files/emit_3_3500_20.ess";
constexpr auto SUCCESS = 0;
constexpr auto FAIL = 1;
} // namespace

int main() {

    launch_GraphBLAS();
    namespace fs = std::experimental::filesystem;

    auto sequences = read_emit_seq(seq_file_name);
    auto is_test_passed = bool(true);

    for (const auto& profile : fs::directory_iterator(chmm_folder)) {

        const auto& path = profile.path();
        auto chmm_name = path.filename().string();

        // Check if file has chmm format and is not "test_chmm.chmm"
        if ((path.extension() == ".chmm") && (chmm_name != "test_chmm.chmm")) {

            auto hmm = read_HMM(path);
            auto non_spec_Viterbi = LA_Viterbi();
            auto spec_Viterbi = LA_Viterbi_spec(path);

            for (const auto& seq : sequences) {
                non_spec_Viterbi.run_Viterbi(hmm, seq);
                spec_Viterbi.run_Viterbi_spec(seq);

                // Check if both results are columns
                auto cols = GrB_Index();
                auto info = GrB_Matrix_ncols(&cols, non_spec_Viterbi.result);
                check_for_error(info);

                if (cols != 1) {
                    std::cerr << "Non spec Viterbi result is not a column!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                info = GrB_Matrix_ncols(&cols, spec_Viterbi.result);
                check_for_error(info);

                if (cols != 1) {
                    std::cerr << "Spec Viterbi result is not a column!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                // Check both results have the same dimension
                auto non_spec_rows = GrB_Index();
                auto spec_rows = GrB_Index();

                info = GrB_Matrix_ncols(&non_spec_rows, non_spec_Viterbi.result);
                check_for_error(info);

                info = GrB_Matrix_ncols(&spec_rows, spec_Viterbi.result);
                check_for_error(info);

                if (non_spec_rows != spec_rows) {
                    std::cerr << "Results' dimensions are incompatible!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                // Check results equality
                auto res_non_spec = HMM::Probability_t();
                auto res_spec = HMM::Probability_t();
                for (size_t i = 0; i < spec_rows; ++i) {
                    info = GrB_Matrix_extractElement_FP32(&res_non_spec, non_spec_Viterbi.result, i,
                                                          0);
                    check_for_error(info);
                    info = GrB_Matrix_extractElement_FP32(&res_spec, spec_Viterbi.result, i, 0);
                    check_for_error(info);

                    if (!HMM::almost_equal(res_non_spec, res_spec) &&
                        (std::isfinite(res_non_spec) || std::isfinite(res_spec))) {
                        std::cerr << "Not equal results at index " << i << ", " << res_non_spec
                                  << " and " << res_spec << ", HMM " << path << '\n';
                        stop_GraphBLAS();
                        return FAIL;
                    }
                }
            }
        }
    }

    stop_GraphBLAS();

    return SUCCESS;
}
