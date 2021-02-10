#include "../LA_Viterbi.h"
#include "../LA_Viterbi_assoc_spec.h"
#include "../LA_Viterbi_spec.h"
#include "../data_reader.h"

#include <experimental/filesystem>
#include <iostream>

namespace {
constexpr auto chmm_folder = "../chmm_files/";
constexpr auto seq_file_name = "../ess_files/emit_3_3500_20.ess";
constexpr auto SUCCESS = 0;
constexpr auto FAIL = 1;

bool is_column(GrB_Matrix res) {
    auto cols = GrB_Index();
    auto info = GrB_Matrix_ncols(&cols, res);
    check_for_error(info);
    return cols == 1;
}

GrB_Index get_rows(GrB_Matrix res) {
    auto rows = GrB_Index(0);
    auto info = GrB_Matrix_nrows(&rows, res);
    check_for_error(info);
    return rows;
}

HMM::Probability_t get_ith_from_column(GrB_Matrix col, GrB_Index i) {
    auto res = HMM::Probability_t();
    auto info = GrB_Matrix_extractElement_FP32(&res, col, i, 0);
    check_for_error(info);
    return res;
}

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
            auto assoc_spec_Viterbi = LA_Viterbi_assoc_spec(path, 2);

            for (const auto& seq : sequences) {
                non_spec_Viterbi.run_Viterbi(hmm, seq);
                spec_Viterbi.run_Viterbi_spec(seq);
                assoc_spec_Viterbi.run_Viterbi_assoc_spec(seq);

                // Check if results are columns

                if (!is_column(non_spec_Viterbi.result)) {
                    std::cerr << "Non spec Viterbi result is not a column!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                if (!is_column(spec_Viterbi.result)) {
                    std::cerr << "Spec Viterbi result is not a column!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                if (!is_column(assoc_spec_Viterbi.result)) {
                    std::cerr << "Associative spec Viterbi result is not a column!\n";
                    stop_GraphBLAS();
                    return FAIL;
                }

                // Check both results have the same dimension
                auto non_spec_rows = get_rows(non_spec_Viterbi.result);
                auto spec_rows = get_rows(spec_Viterbi.result);
                auto assoc_spec_rows = get_rows(assoc_spec_Viterbi.result);

                if (non_spec_rows != spec_rows || spec_rows != assoc_spec_rows) {
                    std::cerr << "Results' dimensions are incompatible!\n"
                              << non_spec_rows << ' ' << spec_rows << ' ' << assoc_spec_rows
                              << '\n';
                    stop_GraphBLAS();
                    return FAIL;
                }

                // Check results equality
                auto res_non_spec = HMM::Probability_t();
                auto res_spec = HMM::Probability_t();
                auto res_assoc_spec = HMM::Probability_t();
                for (size_t i = 0; i < spec_rows; ++i) {
                    res_non_spec = get_ith_from_column(non_spec_Viterbi.result, i);
                    res_spec = get_ith_from_column(spec_Viterbi.result, i);
                    res_assoc_spec = get_ith_from_column(assoc_spec_Viterbi.result, i);

                    if (!HMM::almost_equal(res_non_spec, res_spec) &&
                        (!HMM::almost_equal(res_spec, res_assoc_spec)) &&
                        (std::isfinite(res_non_spec) || std::isfinite(res_spec) ||
                         std::isfinite(res_assoc_spec))) {
                        std::cerr << "Not equal results at index " << i << ", " << res_non_spec
                                  << " vs " << res_spec << " vs " << res_assoc_spec << ", HMM "
                                  << path << '\n';
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
