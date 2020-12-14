#include "data_reader.h"
#include "LA_Viterbi.h"
#include "LA_Viterbi_spec.h"

#include <iostream>
#include <chrono>
#include <experimental/filesystem>
#include <any>


enum class Algorithm_selector { LA, LA_spec };

namespace {
    constexpr std::string_view get_algo_descr(Algorithm_selector alg_sel) {
        if (alg_sel == Algorithm_selector::LA) {
            return "non-specialized version";
        } else if (alg_sel == Algorithm_selector::LA_spec) {
            return "specialized version";
        } else {
            return "Error! Unknown Algorithm_selector";
        }
    }
}


template<int N>
std::chrono::milliseconds benchmark_spec_N_times(
    const std::string& chmm_path,
    const HMM::Seq_vec_t ess)
{
    auto best_time = std::chrono::milliseconds::max();

    auto spec_impl = LA_Viterbi_spec(chmm_path);

    for (size_t i = 0; i < N; ++i) {
        auto cur_iter_time = std::chrono::milliseconds{0};

        for (const auto& seq : ess) {
            spec_impl.run_Viterbi_spec(seq);
        }

        best_time = std::min(best_time, cur_iter_time);
    }

    std::cout << chmm_path << ": best time is " << best_time.count()
        << " msec from " << N << " times\n";
    return best_time;
}


template<int N>
std::chrono::milliseconds benchmark_non_spec_N_times(
    const std::string& chmm_path,
    const HMM::Seq_vec_t ess)
{
    auto best_time = std::chrono::milliseconds::max();

    auto hmm = read_HMM(chmm_path);
    auto non_spec_impl = LA_Viterbi();

    for (size_t i = 0; i < N; ++i) {
        auto cur_iter_time = std::chrono::milliseconds{0};

        for (const auto& seq : ess) {
            non_spec_impl.run_Viterbi(hmm, seq);
        }

        best_time = std::min(best_time, cur_iter_time);
    }

    std::cout << chmm_path << ": best time is " << best_time.count()
        << " msec from " << N << " times\n";
    return best_time;
}


template<int N, Algorithm_selector SEL>
void benchmark_with_chmms_in_folder(
    const std::string& chmm_folder,
    const HMM::Seq_vec_t& ess)
{
    namespace fs = std::experimental::filesystem;
    auto all_time = std::chrono::milliseconds{0};

    std::cout << "\nBenchmarking " << get_algo_descr(SEL)
        << ", chmms from " << chmm_folder << '\n';

    for (const auto& profile : fs::directory_iterator(chmm_folder)) {
        auto path = profile.path();
        auto chmm_name = path.filename().string();
        // Check if file has chmm format and is not "test_chmm.chmm"
        if ((path.extension() == ".chmm") && (chmm_name != "test_chmm.chmm")) {
            if (SEL == Algorithm_selector::LA) {
                all_time += benchmark_non_spec_N_times<N>(path.string(), ess);
            } else if (SEL == Algorithm_selector::LA_spec) {
                all_time += benchmark_spec_N_times<N>(path.string(), ess);
            }
        }
    }

    std::cout << get_algo_descr(SEL) << " best times sum: "
        << all_time.count() << " milliseconds\n\n";
    return;
}
