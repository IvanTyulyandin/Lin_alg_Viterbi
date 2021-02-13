#include "LA_Viterbi.h"
#include "LA_Viterbi_assoc_spec.h"
#include "LA_Viterbi_spec.h"
#include "data_reader.h"

#include <algorithm>
#include <any>
#include <array>
#include <chrono>
#include <experimental/filesystem>
#include <functional>
#include <iostream>

enum class Algo_selector { LA, LA_spec, LA_assoc_spec };

namespace {
constexpr std::string_view get_algo_descr(Algo_selector alg_sel) {
    if (alg_sel == Algo_selector::LA) {
        return "non-specialized version";
    } else if (alg_sel == Algo_selector::LA_spec) {
        return "specialized version";
    } else if (alg_sel == Algo_selector::LA_assoc_spec) {
        return "associative specialized version";
    } else {
        return "Error! Unknown Algorithm_selector";
    }
}
} // namespace

template <size_t N>
std::chrono::milliseconds get_median(const std::function<void(void)>& func,
                                     const std::string& chmm_path) {
    // To prevent out-of-bounds while counting median
    static_assert(N > 1);

    auto results = std::array<std::chrono::milliseconds, N>();

    std::cout << chmm_path << '\n';

    for (size_t i = 0; i < N; ++i) {
        auto iteration_start_time = std::chrono::steady_clock::now();
        func();
        auto cur_time = std::chrono::steady_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - iteration_start_time);
        results[i] = duration;
        std::cout << duration.count() << ' ';
    }

    auto median = std::chrono::milliseconds(0);
    auto mid = results.size() / 2;
    std::sort(results.begin(), results.end());

    if (N % 2) {
        median = results[mid];
    } else {
        median = (results[mid] + results[mid + 1]) / 2;
    }

    std::cout << "\nMedian is " << median.count() << " msec from " << N << " times\n";
    return median;
}

template <int N, Algo_selector SEL, int ASSOC_LEVEL = 1>
void benchmark_with_chmms_in_folder(const std::string& chmm_folder, const HMM::Seq_vec_t& ess) {
    namespace fs = std::experimental::filesystem;
    auto all_time = std::chrono::milliseconds{0};

    std::cout << "\nBenchmarking " << get_algo_descr(SEL) << ", chmms from " << chmm_folder << '\n';

    for (const auto& profile : fs::directory_iterator(chmm_folder)) {
        const auto& path = profile.path();
        auto chmm_name = path.filename().string();
        // Check if file has chmm format and is not "test_chmm.chmm"
        if ((path.extension() == ".chmm") && (chmm_name != "test_chmm.chmm")) {
            if (SEL == Algo_selector::LA) {
                auto hmm = read_HMM(path.string());
                auto non_spec_impl = LA_Viterbi();

                auto non_spec = [&hmm, &ess, &non_spec_impl]() {
                    for (const auto& seq : ess) {
                        non_spec_impl.run_Viterbi(hmm, seq);
                    }
                };
                all_time += get_median<N>(non_spec, chmm_name);

            } else if (SEL == Algo_selector::LA_spec) {
                auto iteration_start_time = std::chrono::steady_clock::now();
                auto spec_impl = LA_Viterbi_spec(path.string());
                auto cur_time = std::chrono::steady_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                    cur_time - iteration_start_time);
                std::cout << "Spec time is " << duration.count() << '\n';
                all_time += duration;

                auto spec = [&ess, &spec_impl]() {
                    for (const auto& seq : ess) {
                        spec_impl.run_Viterbi_spec(seq);
                    }
                };
                all_time += get_median<N>(spec, chmm_name);

            } else if (SEL == Algo_selector::LA_assoc_spec) {
                auto iteration_start_time = std::chrono::steady_clock::now();
                auto spec_assoc_impl = LA_Viterbi_assoc_spec(path.string(), ASSOC_LEVEL);
                auto cur_time = std::chrono::steady_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                    cur_time - iteration_start_time);
                std::cout << "Spec time is " << duration.count() << '\n';
                all_time += duration;

                auto spec_assoc = [&ess, &spec_assoc_impl]() {
                    for (const auto& seq : ess) {
                        spec_assoc_impl.run_Viterbi_assoc_spec(seq);
                    }
                };
                all_time += get_median<N>(spec_assoc, chmm_name);
            }
        }
    }
    if (ASSOC_LEVEL > 1) {
        std::cout << ASSOC_LEVEL << ' ';
    }
    std::cout << get_algo_descr(SEL) << " best times sum: " << all_time.count()
              << " milliseconds\n\n";
    return;
}
