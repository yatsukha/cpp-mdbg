#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <ostream>
#include <string>

#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>

#include <thread_pool_include.hpp>

int main(int argc, char** argv) {
  auto const opts = ::mdbg::command_line_options::parse(argc, argv); 
  auto timer = ::mdbg::timer{};

  if (opts.dry_run) {
    ::std::fprintf(stderr, "### DRY RUN ###\n");
  }

  auto seqs = ::mdbg::load_sequences(opts.input);

  ::std::printf(
    "loaded %lu sequences in %ld ms\n", seqs.size(), timer.reset_ms());

  auto const minimizers = ::mdbg::pick_minimizers(opts.l, opts.d);

  ::std::printf(
    "picked %lu (out of %lu) minimizers (l = %lu, d = %f) in %ld ms\n",
    minimizers.from_hash.size(), 
    static_cast<::std::size_t>(::std::pow(4, opts.l)),
    opts.l, opts.d, timer.reset_ms());
  
  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected(seqs.size());
  ::thread_pool pool{};

  for (::std::size_t i = 0; i < detected.size(); ++i) {
    pool.submit([&seqs, &detected, &minimizers, i]{
      detected[i] = ::mdbg::detect_minimizers(*seqs[i], i, minimizers);
    });
  }

  pool.wait_for_tasks();

  if (!opts.sequences || opts.dry_run) {
    decltype(seqs){}.swap(seqs);
    ::std::printf("cleared sequences, new capacity: %lu\n", seqs.capacity());
  }

  ::std::size_t sum = 0;
  ::std::vector<::std::size_t> stats(detected.size());

  for (::std::size_t i = 0; i < detected.size(); ++i) {
    stats[i] = detected[i].size();
    sum += detected[i].size();
  }

  sum /= detected.size();
  auto const time = timer.reset_ms();
  ::std::sort(stats.begin(), stats.end());

  ::std::printf(
    "detected minimizers in %ld ms\n"
    "minimizers per read:\n"
    "  average:           %lu\n"
    "  median:            %lu\n"
    "  99.9th percentile: %lu\n"
    "  min:               %lu\n",
    time,
    sum,
    stats[stats.size() / 2],
    stats[
      static_cast<::std::size_t>(
        static_cast<double>(stats.size()) * (1.0 - 0.999))],
    stats.front());

  timer.reset_ms();

  auto const& graph =
    ::mdbg::construct(detected.begin(), detected.end(), opts.k);

  ::std::printf(
    "assembled de Bruijn graph (k = %lu) with %lu nodes in %ld ms\n",
    opts.k, graph.size(), timer.reset_ms());

  if (!opts.dry_run) {
    ::std::ofstream out{opts.output};
    ::std::printf("writing...\r"); ::std::fflush(stdout);
    if (!out.is_open()) {
      ::mdbg::terminate("unable to open/create given output file ", opts.output);
    }

    // TODO: problematic usage of sequences
    //       separate tool for loading sequences and graph and outputting gfa?
    ::mdbg::write_gfa(out, graph, seqs, opts);

    ::std::printf(
      "wrote de Bruijn graph to '%s' in %ld ms\n", 
      opts.output.c_str(),
      timer.reset_ms());
  }
}
