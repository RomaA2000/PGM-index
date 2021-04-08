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

#pragma once

#include <sys/stat.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

bool global_verbose = false;

#define IF_VERBOSE(X) if (global_verbose) { X; }
#define OUT_VERBOSE(X) if (global_verbose) { std::cout << "# " << X << std::endl; }

using timer = std::chrono::high_resolution_clock;

template<typename T>
std::string to_metric(const T &x, int digits = 2, bool space = false) {
  static const char *prefix[] = {"y", "z", "a", "f", "p", "n", "µ", "m", "", "k", "M", "G", "T", "P", "E", "Z", "Y"};
  double value = x;
  if (value < 0.)
    value = -value;

  auto log_x = (int) std::log10(value);
  if (log_x > 0)
    log_x = (log_x / 3) * 3;
  else
    log_x = (-log_x + 3) / 3 * (-3);

  value *= std::pow(10, -log_x);
  if (value >= 1000.)
    value /= 1000.0, log_x += 3;

  if (std::fmod(value, 1.0) == 0.0)
    digits = 0;

  constexpr auto prefix_start = -24;
  const char *fmt = space ? (x < 0. ? "-%.*f %s" : "%.*f %s") : (x < 0. ? "-%.*f%s" : "%.*f%s");
  const char *chosen_prefix = prefix[(log_x - prefix_start) / 3];
  std::vector<char> buffer(100);
  int size = std::snprintf(&buffer[0], buffer.size(), fmt, digits, value, chosen_prefix);
  return std::string(&buffer[0], size);
}

template<typename T>
std::vector<T> read_data_binary(const std::string &filename, bool check_sorted) {
  OUT_VERBOSE("Reading " << filename)
  std::vector<T> data;
  try {
    std::fstream in(filename, std::ios::in | std::ios::binary);
    in.exceptions(std::ios::failbit | std::ios::badbit);
    struct stat fs;
    stat(filename.c_str(), &fs);
    size_t file_size = fs.st_size;
    data.resize(file_size / sizeof(T));
    in.read((char *) data.data(), file_size);
  } catch (std::ios_base::failure &e) {
    std::cerr << "Could not read the file. " << e.what() << std::endl;
    exit(1);
  }

  if (size_t(data[0]) == data.size() - 1)
    data.erase(data.begin()); // in some formats, the first element is the size of the dataset, ignore it

  if (check_sorted && !std::is_sorted(data.begin(), data.end())) {
    std::cerr << "Input data must be sorted." << std::endl;
    std::cerr << "Read: [";
    std::copy(data.begin(), std::min(data.end(), data.begin() + 10), std::ostream_iterator<T>(std::cerr, ", "));
    std::cout << "...]." << std::endl;
    exit(1);
  }

  IF_VERBOSE(std::cout << "# Read " << to_metric(data.size()) << " elements: [")
  IF_VERBOSE(std::copy(data.begin(),
                       std::min(data.end() - 1, data.begin() + 5),
                       std::ostream_iterator<T>(std::cout, ", ")))
  IF_VERBOSE(std::cout << "..., " << *(data.end() - 1) << "]" << std::endl)

  return data;
}

template<typename K>
std::vector<char> to_records(const std::vector<K> &keys, size_t record_size) {
  if (record_size < sizeof(K))
    throw std::runtime_error("record_size < key size");
  std::vector<char> out(keys.size() * record_size);
  auto *ptr = out.data();
  for (auto k : keys) {
    *reinterpret_cast<K *>(ptr) = k;
    ptr += record_size;
  }
  return out;
}

#ifdef __GNUG__
#include <memory>
#include <cstdlib>
#include <cxxabi.h>

std::string demangle(const char *name) {
  int status = 0;
  std::unique_ptr<char, void (*)(void *)> res{abi::__cxa_demangle(name, nullptr, nullptr, &status), std::free};
  return status ? name : res.get();
}

#else
std::string demangle(const char* name) { return name; }
#endif

template<typename T>
class RecordIterator {
  const char *ptr;
  size_t step;

 public:
  using value_type = T;
  using difference_type = std::make_signed_t<size_t>;
  using reference = const T &;
  using pointer = const T *;
  using iterator_category = std::random_access_iterator_tag;

  RecordIterator(const char *ptr, size_t step) : ptr(ptr), step(step) {}
  RecordIterator(const RecordIterator &other) : ptr(other.ptr), step(other.step) {}

  RecordIterator &operator=(const RecordIterator &other) {
    ptr = other.ptr;
    step = other.step;
    return *this;
  }

  pointer get() {
    return reinterpret_cast<T *>(ptr);
  }

  reference operator*() const { return *reinterpret_cast<const T *>(ptr); }
  reference operator[](difference_type i) const { return *reinterpret_cast<const T *>(ptr + i * step); }

  RecordIterator &operator++() {
    ptr += step;
    return *this;
  }

  RecordIterator &operator--() {
    ptr -= step;
    return *this;

  }

  RecordIterator operator++(int) {
    RecordIterator it(*this);
    ++(*this);
    return it;
  }

  RecordIterator operator--(int) {
    RecordIterator it(*this);
    --(*this);
    return it;
  }

  RecordIterator operator+(difference_type n) const { return {ptr + n * step, step}; }
  RecordIterator operator-(difference_type n) const { return {ptr - n * step, step}; }

  RecordIterator &operator+=(difference_type n) {
    ptr += n * step;
    return *this;
  }

  RecordIterator &operator-=(difference_type n) {
    ptr -= n * step;
    return *this;
  }

  difference_type operator-(const RecordIterator &it) const {
    assert(step == it.step);
    return (ptr - it.ptr) / it.step;
  }

  bool operator<(const RecordIterator &it) const { return *this - it < 0; }
  bool operator>(const RecordIterator &it) const { return it < *this; }
  bool operator<=(const RecordIterator &it) const { return !(*this > it); }
  bool operator>=(const RecordIterator &it) const { return !(*this < it); }
  bool operator!=(const RecordIterator &it) const { return !(*this == it); }
  bool operator==(const RecordIterator &it) const { return *this - it == 0; }
};

template<typename T>
class NonConstRecordIterator {
  char *ptr;
  size_t step;

 public:
  using value_type = T;
  using difference_type = std::make_signed_t<size_t>;
  using reference = T &;
  using pointer = T *;
  using iterator_category = std::random_access_iterator_tag;

  NonConstRecordIterator(char *ptr, size_t step) : ptr(ptr), step(step) {}
  NonConstRecordIterator(const NonConstRecordIterator &other) : ptr(other.ptr), step(other.step) {}

  NonConstRecordIterator &operator=(const NonConstRecordIterator &other) {
    ptr = other.ptr;
    step = other.step;
    return *this;
  }

  pointer get() {
    return reinterpret_cast<T *>(ptr);
  }

  reference operator*() { return *reinterpret_cast<T *>(ptr); }
  reference operator[](difference_type i) { return *reinterpret_cast<T *>(ptr + i * step); }

  NonConstRecordIterator &operator++() {
    ptr += step;
    return *this;
  }

  NonConstRecordIterator &operator--() {
    ptr -= step;
    return *this;

  }

  NonConstRecordIterator operator++(int) {
    NonConstRecordIterator it(*this);
    ++(*this);
    return it;
  }

  NonConstRecordIterator operator--(int) {
    NonConstRecordIterator it(*this);
    --(*this);
    return it;
  }

  NonConstRecordIterator operator+(difference_type n) const { return {ptr + n * step, step}; }
  NonConstRecordIterator operator-(difference_type n) const { return {ptr - n * step, step}; }

  NonConstRecordIterator &operator+=(difference_type n) {
    ptr += n * step;
    return *this;
  }

  NonConstRecordIterator &operator-=(difference_type n) {
    ptr -= n * step;
    return *this;
  }

  difference_type operator-(const NonConstRecordIterator &it) const {
    assert(step == it.step);
    return (ptr - it.ptr) / it.step;
  }

  bool operator<(const NonConstRecordIterator &it) const { return *this - it < 0; }
  bool operator>(const NonConstRecordIterator &it) const { return it < *this; }
  bool operator<=(const NonConstRecordIterator &it) const { return !(*this > it); }
  bool operator>=(const NonConstRecordIterator &it) const { return !(*this < it); }
  bool operator!=(const NonConstRecordIterator &it) const { return !(*this == it); }
  bool operator==(const NonConstRecordIterator &it) const { return *this - it == 0; }
};

template<typename Class, typename RandomIt>
std::tuple<uint64_t, uint64_t, size_t>
benchmark(RandomIt begin, RandomIt end, const std::vector<typename RandomIt::value_type> &queries) {
  auto t0 = timer::now();
  Class index(begin, end);
  auto t1 = timer::now();
  auto build_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

  auto t2 = timer::now();
  uint64_t cnt = 0;
  for (auto &q : queries) {
    auto range = index.search(q);
    auto lo = begin + range.lo;
    auto hi = begin + range.hi;
    cnt += std::distance(begin, std::lower_bound(lo, hi, q));
  }
  [[maybe_unused]] volatile auto tmp = cnt;
  auto t3 = timer::now();
  auto query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() / queries.size();

  return {build_ms, query_ns, index.size_in_bytes()};
}

template<typename RandomIt>
std::vector<typename RandomIt::value_type>
generate_queries(RandomIt first, RandomIt last, double lookup_ratio, size_t max_queries = 10000000) {
  using value_type = typename RandomIt::value_type;
  auto n = std::distance(first, last);
  auto num_queries = std::min<size_t>(n / 10, max_queries);
  auto num_lookups = size_t(num_queries * lookup_ratio);
  std::uniform_real_distribution<value_type> key_distribution(*first, *std::prev(last));
  std::uniform_real_distribution<value_type> pos_distribution(0, n - 1);
  std::mt19937 generator(std::random_device{}());
  std::vector<value_type> queries;
  queries.reserve(num_queries);

  for (size_t i = 0; i < num_lookups; ++i)
    queries.push_back(first[pos_distribution(generator)]);
  for (size_t i = 0; i < num_queries - num_lookups; ++i)
    queries.push_back(key_distribution(generator));

  std::shuffle(queries.begin(), queries.end(), generator);
  return queries;
}

template<typename T>
struct type_wrapper { using type = T; };

template<typename... Ts, typename TF>
void for_types(TF &&f) { (f(type_wrapper<Ts>{}), ...); }

template<typename K, typename... Args>
size_t benchmark_all(const std::string &filename,
                   std::vector<K> data,
                   double lookup_ratio,
                   const std::string &workload) {
  std::vector<K> queries;
  if (!workload.empty())
    queries = read_data_binary<K>(workload, false);
  else {
    queries = generate_queries(data.begin(), data.end(), lookup_ratio);
    auto m = queries.size();
    OUT_VERBOSE("Generated " << to_metric(m) << " queries, " << to_metric(m * lookup_ratio) << " are lookups")
  }

  size_t min_query_ns = 10000;
  for_types<Args...>([&](auto t) {
    using class_type = typename decltype(t)::type;
    auto name = demangle(typeid(class_type).name());
    auto[build_ms, query_ns, bytes] = benchmark<class_type>(data.begin(), data.end(), queries);
    //std::cout << filename << ",\"" << name << "\", " << build_ms << ", " << bytes << ", " << query_ns << std::endl;
    min_query_ns = std::min(min_query_ns, query_ns);
  });

  return min_query_ns;
}

/*
template<typename K, typename... Args>
void benchmark_all_char(const std::string &filename,
                   const std::vector<char> &data,
                   size_t record_size,
                   double lookup_ratio,
                   const std::string &workload) {
  auto begin = RecordIterator<K>(data.data(), record_size);
  auto end = RecordIterator<K>(data.data() + data.size(), record_size);

  std::vector<K> queries;
  if (!workload.empty())
    queries = read_data_binary<K>(workload, false);
  else {
    queries = generate_queries(begin, end, lookup_ratio);
    auto m = queries.size();
    OUT_VERBOSE("Generated " << to_metric(m) << " queries, " << to_metric(m * lookup_ratio) << " are lookups")
  }

  for_types<Args...>([&](auto t) {
    using class_type = typename decltype(t)::type;
    auto name = demangle(typeid(class_type).name());
    auto[build_ms, query_ns, bytes] = benchmark<class_type>(begin, end, queries);
    std::cout << filename << ",\"" << name << "\", " << build_ms << ", " << bytes << ", " << query_ns << std::endl;
  });
}

template<typename K, typename... Args>
void benchmark_presaved(size_t record_size, double lookup_ratio) {
  std::ifstream input_file_data("");
  std::vector<K> data;
  K value2;
  data.push_back(1.0001);
  for (size_t i = 1; i < 1000000; ++i) {
    data.push_back(data[i-1]*data[0]);
  }
  K value1;
  std::vector<K> queries;
  queries = generate_queries(data.begin(), data.end(), lookup_ratio);
  auto m = queries.size();
  OUT_VERBOSE("Generated " << to_metric(m) << " queries, " << to_metric(m * lookup_ratio) << " are lookups")


  for_types<Args...>([&](auto t) {
    using class_type = typename decltype(t)::type;
    auto name = demangle(typeid(class_type).name());
    auto[build_ms, query_ns, bytes] = benchmark<class_type>(data.begin(), data.end(), queries);
    std::cout << name << "\"," << build_ms << "," << bytes << "," << query_ns << std::endl;
  });
}
*/


/*
template<typename K>
void benchmark_binary_search_char(const std::string &filename,
                             const std::vector<char> &data,
                             size_t record_size,
                             double lookup_ratio) {
  auto begin = RecordIterator<K>(data.data(), record_size);
  auto end = RecordIterator<K>(data.data() + data.size(), record_size);
  auto queries = generate_queries(begin, end, lookup_ratio);
  auto query_ns = benchmark_simple_binary(begin, end, queries);
  auto t_start = timer::now();
  std::vector<char> data_optimized = data;

  auto begin_optimized = NonConstRecordIterator<K>(data_optimized.data(), record_size);
  auto end_optimized = NonConstRecordIterator<K>(data_optimized.data() + data_optimized.size(), record_size);

  eytzinger_maker<NonConstRecordIterator<K>, RecordIterator<K>, K>(begin_optimized, end_optimized, begin);
  auto t_finish = timer::now();
  auto building = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();

  queries = generate_queries(begin, end, lookup_ratio);
  auto query_ns_optimized = benchmark_binary(begin_optimized, end_optimized, queries);
  std::cout << filename << ",\"" << "Simple_binary_search" << "\", " << 0 << ", " << query_ns << std::endl;
  std::cout << filename << ",\"" << "Optimized_binary_search" << "\", " << building << ", " << query_ns_optimized << std::endl;
}
*/


template<typename RandomIt1, typename RandomIt2>
size_t
eytzinger_maker(RandomIt1 begin, RandomIt1 end, RandomIt2 begin_prev_data, size_t i = 0, size_t k = 1) {
  using K = typename RandomIt1::value_type;
  if (k <= end - begin) {
    i = eytzinger_maker<RandomIt1, RandomIt2>(begin, end, begin_prev_data, i, 2 * k);
    begin[k - 1] = K(begin_prev_data[i++]);
    i = eytzinger_maker<RandomIt1, RandomIt2>(begin, end, begin_prev_data, i, 2 * k + 1);
  }
  return i;
}

template<typename RandomIt>
RandomIt
eytzinger(RandomIt begin, RandomIt end, typename RandomIt::value_type q) {
  size_t k = 1u;
  size_t n = end - begin;
  while (k <= n) {
    __builtin_prefetch(static_cast<void *>(&(begin[0]) + k * 16));
    k = 2 * k + (begin[k - 1] < q);
  }
  k >>= __builtin_ffs(~k);
  return begin + (k - 1);
}


template<typename RandomIt>
uint64_t
benchmark_optimized_binary(RandomIt begin, RandomIt end, const std::vector<typename RandomIt::value_type> &queries) {
  auto t2 = timer::now();
  uint64_t cnt = 0;
  for (auto &q : queries) {
    cnt += std::distance(begin, eytzinger(begin, end, q));
  }
  [[maybe_unused]] volatile auto tmp = cnt;
  auto t3 = timer::now();
  auto query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() / queries.size();
  return query_ns;
}

template<typename RandomIt>
uint64_t
benchmark_simple_binary(RandomIt begin, RandomIt end, const std::vector<typename RandomIt::value_type> &queries) {
  auto t2 = timer::now();
  uint64_t cnt = 0;
  for (auto &q : queries) {
    cnt += std::distance(begin, std::lower_bound(begin, end, q));
  }
  [[maybe_unused]] volatile auto tmp = cnt;
  auto t3 = timer::now();
  auto query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() / queries.size();
  return query_ns;
}


template<typename K>
size_t benchmark_binary(const std::string &filename,
                                  const std::vector<K> &data,
                                  double lookup_ratio) {
  std::vector<K> queries = generate_queries(data.begin(), data.end(), lookup_ratio);

  auto query_ns = benchmark_simple_binary(data.begin(), data.end(), queries);
  //std::cout << filename << ",\"" << "Simple_binary_search" << "\", " << 0 << ", " << query_ns << std::endl;
  
  auto t_start = timer::now();
  std::vector<K> data_optimized = data;
  eytzinger_maker(data_optimized.begin(), data_optimized.end(), data.begin());
  auto t_finish = timer::now();
  auto building = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();

  auto query_ns_optimized = benchmark_optimized_binary(data_optimized.begin(), data_optimized.end(), queries);
  //std::cout << filename << ",\"" << "Optimized_binary_search" << "\", " << building << ", " << query_ns_optimized << std::endl;

  return query_ns_optimized;
}