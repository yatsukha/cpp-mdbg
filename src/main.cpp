#include "mdbg/simplification.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>

#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>
#include <mdbg/trio_binning.hpp>

#include <thread_pool_include.hpp>

void construct_from(
  ::mdbg::sequences_t const& seqs,
  ::thread_pool& pool,
  ::mdbg::command_line_options const& opts,
  char const* print_prefix = ""
) noexcept {
  auto timer = ::mdbg::timer{};

  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected(seqs.size());
  ::std::vector<::std::future<bool>> detected_futures;

  for (::std::size_t i = 0; i < detected.size(); ++i) {
    detected_futures.emplace_back(
      pool.submit([&seqs, &detected, &opts, i]{
        detected[i] = ::mdbg::detect_minimizers(*seqs[i], i, opts);
      })
    );
  }

  // pool.wait_for_tasks();

  ::std::size_t sum = 0;
  ::std::vector<::std::size_t> stats(detected.size());

  for (::std::size_t i = 0; i < detected.size(); ++i) {
    detected_futures[i].get();
    stats[i] = detected[i].size();
    sum += detected[i].size();
  }

  sum /= detected.size();
  auto const time = timer.reset_ms();
  ::std::sort(stats.begin(), stats.end());

  ::std::printf(
    "%sdetected minimizers in %ld ms\n"
    "%sminimizers per read:\n"
    "  average:           %lu\n"
    "  median:            %lu\n"
    "  99.9th percentile: %lu\n"
    "  min:               %lu\n",
    print_prefix, time,
    print_prefix,
    sum,
    stats[stats.size() / 2],
    stats[
      static_cast<::std::size_t>(
        static_cast<double>(stats.size()) * (1.0 - 0.999))],
    stats.front());
  ::std::fflush(stdout);
  
  if (opts.analysis) {
    return;
  }

  timer.reset_ms();

  auto const& graph =
    ::mdbg::construct(detected.begin(), detected.end(), opts);

  ::std::printf(
    "%sassembled de Bruijn graph (k = %lu) with %lu nodes in %ld ms\n",
    print_prefix, opts.k, graph.size(), timer.reset_ms());
  ::std::fflush(stdout);

  ::mdbg::simplified_graph_t simplified;

  if (opts.unitigs) {
    simplified = ::mdbg::simplify(graph);
    ::std::printf(
      "%ssimplified to %lu nodes in %ld ms\n",
      print_prefix, simplified.size(), timer.reset_ms());
    ::std::fflush(stdout);
  }

  if (!opts.dry_run) {
    ::std::ofstream out{opts.output_prefix};
    if (!out.is_open()) {
      ::mdbg::terminate("unable to open/create given output file ", opts.output_prefix);
    }

    if (opts.unitigs) {
      ::mdbg::write_gfa(out, simplified, seqs, opts);
    } else {
      ::mdbg::write_gfa(out, graph, seqs, opts);
    }

    ::std::printf(
      "%swrote de Bruijn graph to '%s' in %ld ms\n",
      print_prefix,
      opts.output_prefix.c_str(),
      timer.reset_ms());
    ::std::fflush(stdout);
  }
}


int main(int argc, char** argv) {
  auto opts = ::mdbg::command_line_options::parse(argc, argv); 
  auto timer = ::mdbg::timer{};

  if (opts.dry_run) {
    ::std::fprintf(stderr, "### DRY RUN  ###\n");
  }

  if (opts.analysis) {
    ::std::fprintf(stderr, "### ANALYSIS ###\n");
  }

  auto seqs = ::mdbg::load_sequences(opts.input);

  ::std::printf(
    "loaded %lu sequences in %ld ms\n", seqs.size(), timer.reset_ms());

  ::thread_pool pool{opts.threads};

  if (!opts.trio_binning) {
    opts.output_prefix += ".gfa";
    ::construct_from(seqs, pool, opts);
    ::std::exit(EXIT_SUCCESS);
  }

  ::std::vector<::mdbg::kmer_counts_t> counts(2);

  pool.submit([&counts, &opts]{
    counts[0] = ::mdbg::count_kmers(
      ::mdbg::load_sequences(opts.trio_binning->input_0),
      opts);
  });

  pool.submit([&counts, &opts]{
    counts[1] = ::mdbg::count_kmers(
      ::mdbg::load_sequences(opts.trio_binning->input_1),
      opts);
  });

  pool.wait_for_tasks();

  ::std::printf(
    "loaded and counted kmers of haplotypes in %ld ms\n", timer.reset_ms());

  ::mdbg::reduce_to_unique_kmers(counts[0], counts[1]);
  
  ::std::printf(
    "reduced to unique kmers (%lu, %lu) in %ld ms\n",
    counts[0].size(), counts[1].size(), timer.reset_ms());

  // TODO: speed this part up with threadpool
  auto const& filtered = ::mdbg::filter_reads(pool, seqs, counts, opts);

  ::std::printf(
    "binned reads (0 - %lu, 1 - %lu) in %ld ms\n", 
    filtered.first.size(), filtered.second.size(), timer.reset_ms());

  auto opts_0 = opts;
  opts_0.output_prefix += ".0.gfa";

  pool.submit([&filtered, &pool, &opts_0]{
    ::construct_from(filtered.first, pool, opts_0, "[haplotype 0] ");
  });

  auto opts_1 = opts;
  opts_1.output_prefix += ".1.gfa";

  pool.submit([&filtered, &pool, &opts_1]{
    ::construct_from(filtered.second, pool, opts_1, "[haplotype 1] ");
  });

  pool.wait_for_tasks();
}
