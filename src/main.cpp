#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>

#include <mdbg/opt.hpp>
#include <mdbg/io/parser.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/graph/simplification.hpp>
#include <mdbg/graph/construction.hpp>
#include <mdbg/trio_binning/trio_binning.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/task_group.h>

int main(int argc, char** argv) {
  auto opts = ::mdbg::command_line_options::parse(argc, argv);
  auto timer = ::mdbg::timer{};

  if (opts.threads) {
    ::tbb::global_control max_parallelism{
      ::tbb::global_control::max_allowed_parallelism, opts.threads};
  }

  if (opts.dry_run) {
    ::std::fprintf(stderr, "### DRY RUN  ###\n");
  }

  if (opts.analysis) {
    ::std::fprintf(stderr, "### ANALYSIS ###\n");
  }
  
  using processed_pair_t = ::std::pair<::std::string, ::mdbg::read_minimizers_t>;

  ::std::vector<::std::unique_ptr<processed_pair_t>> processed;
  ::mdbg::graph::de_bruijn_graph_t graph;
  ::tbb::task_group tg;
  ::std::mutex tg_lock;

  ::mdbg::io::fasta_constumer consumer = 
    [&processed, &graph, &tg, &tg_lock, &opts](auto&&, auto&& seq) {
      processed.emplace_back(
        ::std::make_unique<processed_pair_t>(
          processed_pair_t{::std::forward<::std::string>(seq), {}}));
      auto const index = processed.size() - 1;

      tg_lock.lock();
      tg.run([&tg, &tg_lock, ptr = processed.back().get(), &graph, &opts, index]{
        ptr->second = ::mdbg::detect_minimizers(ptr->first, index, opts);

        if (!opts.analysis) {
          tg_lock.lock();
          tg.run([&graph, ptr, &opts]{
            ::mdbg::graph::construct(graph, ptr->second, opts);
          });
          tg_lock.unlock();
        }
      });
      tg_lock.unlock();
    };


  ::mdbg::io::parse_fasta(opts.input.c_str(), consumer);
  tg.wait();

  ::std::printf(
    "loaded, processed and assembled a graph from "
    "sequences from %lu sequences in %ld ms\n", 
    processed.size(), timer.reset_ms());


  opts.output_prefix += ".gfa";

  ::std::vector<::std::size_t> stats(processed.size());

  for (::std::size_t i = 0; i < stats.size(); ++i) {
    stats[i] = processed[i]->second.size();
  }

  auto const time = timer.reset_ms();
  ::std::sort(stats.begin(), stats.end());
  auto const get_percentile = [&stats](double const percentile) {
    return stats[static_cast<::std::size_t>(
      static_cast<double>(stats.size()) * (1.0 - percentile))];
  };

  ::std::printf(
    "calculated statistics in %ld ms\n"
    "minimizers per read:\n"
    "  median:             %lu\n"
    "  90th    percentile: %lu\n"
    "  99th    percentile: %lu\n"
    "  99.9th  percentile: %lu\n"
    "  99.99th percentile: %lu\n"
    "  min:                %lu\n",
    timer.reset_ms(),
    get_percentile(0.5),
    get_percentile(0.9),
    get_percentile(0.99),
    get_percentile(0.999),
    get_percentile(0.9999),
    stats.front());
  ::std::fflush(stdout);
  
  if (opts.analysis) {
    ::std::exit(EXIT_SUCCESS);
  }

  ::std::printf(
    "assembled de Bruijn graph (k = %lu) with %lu node(s)\n",
    opts.k, graph.size());
  ::std::fflush(stdout);

  auto const simplified = ::mdbg::graph::simplify(graph);

  ::std::printf(
    "simplified to %lu node(s) in %ld ms\n",
    simplified.size(), timer.reset_ms());
  ::std::fflush(stdout);

  if (!opts.dry_run) {
    ::std::ofstream out{opts.output_prefix};
    if (!out.is_open()) {
      ::mdbg::terminate("unable to open/create given output file ", opts.output_prefix);
    }

    ::mdbg::graph::write_gfa(
      out, 
      simplified, 
      [&processed](auto&& i) -> std::string const& { return processed[i]->first; },
      opts);
    
    ::std::printf(
      "wrote de Bruijn graph to '%s' in %ld ms\n",
      opts.output_prefix.c_str(),
      timer.reset_ms());
    ::std::fflush(stdout);
  }
}
