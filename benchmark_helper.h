#include "LA_Viterbi.h"
#include "LA_Viterbi_assoc_spec.h"
#include "LA_Viterbi_spec.h"
#include "data_reader.h"

#include <any>
#include <chrono>
#include <experimental/filesystem>
#include <iostream>

enum class Algorithm_selector { LA, LA_spec, LA_assoc_spec };
constexpr size_t ASSOC_LEVEL = 2;

namespace {
constexpr std::string_view get_algo_descr(Algorithm_selector alg_sel) {
    if (alg_sel == Algorithm_selector::LA) {
        return "non-specialized version";
    } else if (alg_sel == Algorithm_selector::LA_spec) {
        return "specialized version";
    } else if (alg_sel == Algorithm_selector::LA_assoc_spec) {
        return "associative specialized version";
    } else {
        return "Error! Unknown Algorithm_selector";
    }
}
} // namespace

template <int N>
std::chrono::milliseconds benchmark_non_spec_N_times(const std::string& chmm_path,
                                                     const HMM::Seq_vec_t& ess) {
    auto best_time = std::chrono::milliseconds::max();

    auto hmm = read_HMM(chmm_path);
    auto non_spec_impl = LA_Viterbi();

    for (size_t i = 0; i < N; ++i) {
        auto iteration_start_time = std::chrono::steady_clock::now();

        for (const auto& seq : ess) {
            non_spec_impl.run_Viterbi(hmm, seq);
        }

        auto cur_time = std::chrono::steady_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - iteration_start_time);
        best_time = std::min(best_time, duration);
    }

    std::cout << chmm_path << ": best time is " << best_time.count() << " msec from " << N
              << " times\n";
    return best_time;
}

template <int N>
std::chrono::milliseconds benchmark_spec_N_times(const std::string& chmm_path,
                                                 const HMM::Seq_vec_t& ess) {
    auto best_time = std::chrono::milliseconds::max();

    auto spec_impl = LA_Viterbi_spec(chmm_path);

    for (size_t i = 0; i < N; ++i) {
        auto iteration_start_time = std::chrono::steady_clock::now();

        for (const auto& seq : ess) {
            spec_impl.run_Viterbi_spec(seq);
        }

        auto cur_time = std::chrono::steady_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - iteration_start_time);
        best_time = std::min(best_time, duration);
    }

    std::cout << chmm_path << ": best time is " << best_time.count() << " msec from " << N
              << " times\n";
    return best_time;
}

template <int N>
std::chrono::milliseconds benchmark_assoc_spec_N_times(const std::string& chmm_path,
                                                       const HMM::Seq_vec_t& ess,
                                                       const size_t level) {
    auto best_time = std::chrono::milliseconds::max();

    auto spec_impl = LA_Viterbi_assoc_spec(chmm_path, level);

    for (size_t i = 0; i < N; ++i) {
        auto iteration_start_time = std::chrono::steady_clock::now();

        for (const auto& seq : ess) {
            spec_impl.run_Viterbi_assoc_spec(seq);
        }

        auto cur_time = std::chrono::steady_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - iteration_start_time);
        best_time = std::min(best_time, duration);
    }

    std::cout << chmm_path << ": best time is " << best_time.count() << " msec from " << N
              << " times\n";
    return best_time;
}

template <int N, Algorithm_selector SEL>
void benchmark_with_chmms_in_folder(const std::string& chmm_folder, const HMM::Seq_vec_t& ess) {
    namespace fs = std::experimental::filesystem;
    auto all_time = std::chrono::milliseconds{0};

    std::cout << "\nBenchmarking " << get_algo_descr(SEL) << ", chmms from " << chmm_folder << '\n';

    for (const auto& profile : fs::directory_iterator(chmm_folder)) {
        const auto& path = profile.path();
        auto chmm_name = path.filename().string();
        // Check if file has chmm format and is not "test_chmm.chmm"
        if ((path.extension() == ".chmm") && (chmm_name != "test_chmm.chmm")) {
            if (SEL == Algorithm_selector::LA) {
                all_time += benchmark_non_spec_N_times<N>(path.string(), ess);
            } else if (SEL == Algorithm_selector::LA_spec) {
                all_time += benchmark_spec_N_times<N>(path.string(), ess);
            } else if (SEL == Algorithm_selector::LA_assoc_spec) {
                all_time += benchmark_assoc_spec_N_times<N>(path.string(), ess, ASSOC_LEVEL);
            }
        }
    }

    std::cout << get_algo_descr(SEL) << " best times sum: " << all_time.count()
              << " milliseconds\n\n";
    return;
}
