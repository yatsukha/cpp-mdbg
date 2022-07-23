#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>

#include <mdbg/simplification.hpp>
#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>
#include <mdbg/trio_binning.hpp>

#include <thread_pool_include.hpp>

#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>

void construct_from(
  ::mdbg::sequences_t const& seqs,
  ::mdbg::command_line_options const& opts,
  char const* print_prefix = ""
) noexcept {
  auto timer = ::mdbg::timer{};

  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected(seqs.size());
  ::std::vector<::std::size_t> stats(detected.size());

  ::tbb::parallel_for(
    ::tbb::blocked_range<::std::size_t>(0, seqs.size()),
    [&detected, &seqs, &opts, &stats](
      ::tbb::blocked_range<::std::size_t> const& r
    ) {
      for (auto i = r.begin(); i != r.end(); ++i) {
        detected[i] = ::mdbg::detect_minimizers(*seqs[i], i, opts);
        stats[i] = detected[i].size();
      }
    });

  auto const time = timer.reset_ms();
  ::std::sort(stats.begin(), stats.end());
  auto const get_percentile = [&stats](double const percentile) {
    return stats[static_cast<::std::size_t>(
      static_cast<double>(stats.size()) * (1.0 - percentile))];
  };

  ::std::printf(
    "%sdetected minimizers in %ld ms\n"
    "%sminimizers per read:\n"
    "  median:             %lu\n"
    "  90th    percentile: %lu\n"
    "  99th    percentile: %lu\n"
    "  99.9th  percentile: %lu\n"
    "  99.99th percentile: %lu\n"
    "  min:                %lu\n",
    print_prefix, time,
    print_prefix,
    get_percentile(0.5),
    get_percentile(0.9),
    get_percentile(0.99),
    get_percentile(0.999),
    get_percentile(0.9999),
    stats.front());
  ::std::fflush(stdout);
  
  if (opts.analysis) {
    return;
  }

  timer.reset_ms();

  ::mdbg::de_bruijn_graph_t graph;

  ::tbb::parallel_for(
    ::tbb::blocked_range<::std::size_t>(0, detected.size()), 
    [&detected, &graph, &opts](::tbb::blocked_range<::std::size_t> const& r) {
      for (auto i = r.begin(); i != r.end(); ++i) {
        ::mdbg::construct(graph, detected[i], opts);
      }
    });

  ::std::printf(
    "%sassembled de Bruijn graph (k = %lu) with %lu node(s) in %ld ms\n",
    print_prefix, opts.k, graph.size(), timer.reset_ms());
  ::std::fflush(stdout);

  ::mdbg::simplified_graph_t simplified;

  if (opts.unitigs) {
    simplified = ::mdbg::simplify(graph);
    ::std::printf(
      "%ssimplified to %lu node(s) in %ld ms\n",
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
      ::std::fprintf(
        stderr, 
        "non-simplified graph output is deprecated and may not work as expected; "
        "use -u\n");
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
  ::tbb::global_control max_parallelism{
    ::tbb::global_control::max_allowed_parallelism, pool.get_thread_count()};

  if (!opts.trio_binning) {
    opts.output_prefix += ".gfa";
    ::construct_from(seqs, opts);
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

  pool.submit([&filtered, &opts_0]{
    ::construct_from(filtered.first, opts_0, "[haplotype 0] ");
  });

  auto opts_1 = opts;
  opts_1.output_prefix += ".1.gfa";

  pool.submit([&filtered, &opts_1]{
    ::construct_from(filtered.second, opts_1, "[haplotype 1] ");
  });

  pool.wait_for_tasks();
}
