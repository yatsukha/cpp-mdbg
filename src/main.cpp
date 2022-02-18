#include "mdbg/io.hpp"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

int main(int argc, char** argv) {
  auto const opts = ::mdbg::command_line_options::parse(argc, argv); 
  auto timer = ::mdbg::timer{};

  auto const seqs = ::mdbg::load_sequences(opts.input);

  ::std::printf(
    "loaded %lu sequences in %ld ms\n", seqs.size(), timer.reset_ms());

  auto const minimizers = ::mdbg::pick_minimizers(opts.l, opts.d);

  ::std::printf(
    "picked %lu minimizers of length %lu in %ld ms\n",
    minimizers.from_hash.size(), opts.l, timer.reset_ms());

  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected(seqs.size());

  ::std::size_t sum = 0;
  ::std::vector<::std::size_t> statistics(seqs.size());

  for (::std::size_t i = 0; i < seqs.size(); ++i) {
    detected[i] = ::mdbg::detect_minimizers(*seqs[i], i, minimizers);
    statistics[i] = detected[i].size();
    sum += detected[i].size();
  }

  sum /= seqs.size();
  auto const time = timer.reset_ms();
  ::std::sort(statistics.begin(), statistics.end());

  ::std::printf(
    "detected minimizers in %ld ms\n"
    "minimizers per read (l = %lu, d = %f):\n"
    "  average:         %lu\n"
    "  median:          %lu\n"
    "  99th percentile: %lu\n",
    time,
    opts.l, opts.d,
    sum,
    statistics[statistics.size() / 2],
    statistics[
      static_cast<::std::size_t>(
        static_cast<double>(statistics.size()) * (1.0 - 0.99))]
  );

  timer.reset_ms();
  auto graph = ::mdbg::construct(detected, opts.k);

  ::std::printf(
    "assembled de Bruijn graph (k = %lu) with %lu nodes in %ld ms\n",
    opts.k, graph.size(), timer.reset_ms());
}
