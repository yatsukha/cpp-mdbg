#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>

#include <fast_include.hpp>
#include <biosoup_include.hpp>

#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/fast_custom_factory.hpp>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

// TODO: program options - https://github.com/jarro2783/cxxopts/tree/v3.0.0
int main(int argc, char** argv) {
  if (argc != 2) {
    ::mdbg::terminate("Required arguments: <input reads>");
  }

  if (!::std::filesystem::exists(argv[1])) {
    ::mdbg::terminate("Could not locate given file: ", argv[1]);
  }

  auto timer = ::mdbg::timer{};

  // load the sequences, ty @tbrekalo
  auto parser =
      fast::CreateFastaParser<::std::string>(argv[1]);
  auto const seqs = parser.Parse();

  ::std::printf(
    "loaded %lu sequences in %ld ms\n", seqs.size(), timer.reset_ms());

  auto constexpr l = 7;
  auto constexpr d = 0.008;

  auto const minimizers = ::mdbg::pick_minimizers(l, d);

  ::std::printf(
    "picked %lu minimizers of length %d in %ld ms\n",
    minimizers.from_hash.size(), l, timer.reset_ms());

  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected;
  detected.reserve(seqs.size());

  ::std::size_t sum = 0;
  ::std::vector<::std::size_t> statistics(seqs.size());
  for (::std::size_t i = 0; i < seqs.size(); ++i) {
    detected.emplace_back(::mdbg::detect_minimizers(*seqs[i], i, minimizers));
    statistics[i] = detected.back().size();
    sum += detected.back().size();
  }

  sum /= seqs.size();
  auto const time = timer.reset_ms();
  ::std::sort(statistics.begin(), statistics.end());

  ::std::printf(
    "detected minimizers in %ld ms\n"
    "minimizers per read (l = %d, d = %f):\n"
    "  average:         %lu\n"
    "  median:          %lu\n"
    "  99th percentile: %lu\n",
    time,
    l, d,
    sum,
    statistics[statistics.size() / 2],
    statistics[
      static_cast<::std::size_t>(
        static_cast<double>(statistics.size()) * (1.0 - 0.99))
    ]
  );

  timer.reset_ms();

  auto constexpr k = 5;
  auto graph = ::mdbg::construct(detected, k);

  ::std::printf(
    "assembled de Bruijn graph (k = %d) with %lu nodes in %ld ms\n",
    k, graph.size(), timer.reset_ms());
}
