// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2021 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "benchmark.hpp"
#include "args.hxx"
#include "pgm/pgm_index.hpp"
#include "pgm/pgm_index_variants.hpp"

#include <fstream>
#include <functional>
#include <utility>

#define FOR_EACH_EPSILON(C, K) C<K, 8>, C<K, 16>, C<K, 32>, C<K, 64>, C<K, 128>, C<K, 256>, \
                               C<K, 512>, C<K, 1024>, C<K, 2048>, C<K, 4096>
#define FOR_EACH_EPSILON_FLOATING(C, K) C<K, 8, 4, K>, C<K, 16, 4, K>, C<K, 32, 4, K>, C<K, 64, 4, K>, C<K, 128, 4, K>, C<K, 256, 4, K>, \
                               C<K, 512, 4, K>, C<K, 1024, 4, K>, C<K, 2048, 4, K>, C<K, 4096, 4, K>

#define PGM_CLASSES(K) FOR_EACH_EPSILON(pgm::PGMIndex, K)
#define BPGM_CLASSES(K) FOR_EACH_EPSILON(pgm::BucketingPGMIndex, K)
#define EFPGM_CLASSES(K) FOR_EACH_EPSILON(pgm::EliasFanoPGMIndex, K)
#define CPGM_CLASSES(K) FOR_EACH_EPSILON(pgm::CompressedPGMIndex, K)

#define PGM_CLASSES_FLOATING(K) FOR_EACH_EPSILON_FLOATING(pgm::PGMIndex, K)
#define BPGM_CLASSES_FLOATING(K) FOR_EACH_EPSILON(pgm::BucketingPGMIndex, K)
#define EFPGM_CLASSES_FLOATING(K) FOR_EACH_EPSILON(pgm::EliasFanoPGMIndex, K)
#define CPGM_CLASSES_FLOATING(K) FOR_EACH_EPSILON(pgm::CompressedPGMIndex, K)

#define ALL_CLASSES(K) PGM_CLASSES(K), BPGM_CLASSES(K), EFPGM_CLASSES(K), CPGM_CLASSES(K)

#define ALL_CLASSES_FLOATING(K) PGM_CLASSES_FLOATING(K)
//BPGM_CLASSES_FLOATING(K), EFPGM_CLASSES_FLOATING(K), CPGM_CLASSES_FLOATING(K)

template<typename K>
void read_ints_helper(args::PositionalList<std::string> &files,
                      size_t record_size,
                      int lookup_ratio,
                      const std::string &workload) {
    OUT_VERBOSE("Running with " << sizeof(K) << "-byte keys + " << record_size - sizeof(K) << "-byte values")
    for (const auto &file : files.Get()) {
        auto data = to_records(read_data_binary<K>(file, true), record_size);
        auto filename = file.substr(file.find_last_of("/\\") + 1);
        benchmark_all<K, ALL_CLASSES(K)>(filename, data, record_size, lookup_ratio, workload);
    }
}


int main(int argc, char **argv) {
    using namespace args;
    ArgumentParser p("Benchmark for the PGM-index library.");
    p.helpParams.flagindent = 2;
    p.helpParams.helpindent = 25;
    p.helpParams.progindent = 0;
    p.helpParams.descriptionindent = 0;

    CompletionFlag completion(p, {"complete"});
    HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
    Flag verbose(p, "", "Verbose output", {'v', "verbose"});
    ValueFlag<size_t> value_size(p, "bytes", "Size of the values associated to keys", {'V', "values"}, 0);

    Group g1(p, "QUERY WORKLOAD OPTIONS (mutually exclusive):", Group::Validators::AtMostOne);
    ValueFlag<int> ratio(g1, "ratio", "Random workload with the given lookup ratio", {'r', "ratio"}, 0.333);
    ValueFlag<std::string> workload(g1, "file", "Custom workload file. Obeys the format of input files", {'w', "workload"});

    Group g2(p, "INPUT DATA OPTIONS (mutually exclusive):", Group::Validators::Xor, Options::Required);
    ValueFlag<size_t> synthetic(g2, "size", "Generate synthetic data of the given size", {'s', "synthetic"}, 1000000000);
    Flag u64(g2, "", "Input files contain unsigned 64-bit ints", {'U', "u64"});
    Flag i64(g2, "", "Input files contain signed 64-bit ints", {'I', "i64"});
    PositionalList<std::string> files(p, "file", "The input files");

    try {
        p.ParseCLI(argc, argv);
    }
    catch (args::Completion &e) {
        std::cout << e.what();
        return 0;
    }
    catch (args::Help &) {
        std::cout << p;
        return 0;
    }
    catch (args::Error &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << p;
        return 1;
    }

    if (ratio.Get() < 0 || ratio.Get() > 1) {
        std::cerr << "Argument to --" << ratio.GetMatcher().GetLongOrAny().str() << " must be between 0.0 and 1.0.";
        return 1;
    }

    if (synthetic && synthetic.Get() < 1000) {
        std::cerr << "Argument to --" << synthetic.GetMatcher().GetLongOrAny().str() << " must be greater than 1000.";
        return 1;
    }

    global_verbose = verbose.Get();
    std::cout << "dataset,class_name,build_ms,bytes,query_ns" << std::endl;

    if (synthetic) {
        auto record_size = value_size.Get() + sizeof(uint64_t);
        auto n = synthetic.Get();
        std::mt19937 generator(std::random_device{}());
        auto gen = [&](auto distribution) {
            std::vector<long double> out(n);
            std::generate(out.begin(), out.end(), [&] { return distribution(generator); });
            std::sort(out.begin(), out.end());
            return out;
        };
        std::vector<std::pair<std::string, std::function<std::vector<long double>()>>> distributions = {
            {"exponential", [n]() {
                long double base = 1.00001;
                std::vector<long double> data;
                data.push_back(base);
                for (int i = 1; i < n; i++)
                {
                    data.push_back(data[i-1] * base);
                }
                static int count = 0;
                count++;
                if (count == 1)
                {
                    std::cout << "Maximum of exponent: " << data[n-1] << std::endl;
                }
                return data;} },
/*
            {"uniform_dense", std::bind(gen, std::uniform_real_distribution<long double>(0, n * 1000))},
            {"uniform_sparse", std::bind(gen, std::uniform_real_distribution<long double>(0, n * n))},
            {"binomial", std::bind(gen, std::binomial_distribution<int>(1ull << 50))},
            {"negative_binomial", std::bind(gen, std::negative_binomial_distribution<int>(1ull << 50, 0.3))},
            {"geometric", std::bind(gen, std::geometric_distribution<int>(1e-10))},
*/
        };
        OUT_VERBOSE("Generating " << to_metric(n) << " elements (8-byte keys + " << value_size.Get() << "-byte values)")
        
        double mean_advantage = 0.0;
        size_t num_tests = 5;

        for (size_t i = 0; i < num_tests; i++)
        {
            for (auto&[name, gen_data] : distributions) {
              double min_pgm_ns = benchmark_all<long double, ALL_CLASSES_FLOATING(long double)>(name, gen_data(), ratio.Get(), workload.Get());
              double min_bin_ns = benchmark_binary<long double>(name, gen_data(), ratio.Get());

              mean_advantage += min_pgm_ns / min_bin_ns;
            }
        }

        mean_advantage *= 100;
        mean_advantage /= num_tests;
        std::cout << "mean_advantage = " << mean_advantage << "%" << std::endl;
    }


//    if (i64.Get())
//        read_ints_helper<int64_t>(files, value_size.Get() + sizeof(int64_t), ratio.Get(), workload.Get());
//    if (u64.Get())
//        read_ints_helper<uint64_t>(files, value_size.Get() + sizeof(uint64_t), ratio.Get(), workload.Get());

    return 0;
}