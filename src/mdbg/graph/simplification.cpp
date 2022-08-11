#include <cstdio>
#include <mutex>

#include <mdbg/graph/simplification.hpp>
#include <mdbg/graph/construction.hpp>

#include <tbb/task_group.h>

namespace mdbg::graph {


  simplified_graph_t::mapped_type unitig(
    de_bruijn_graph_t const& dbg,
    detail::compact_minimizer const& starting_minimizer
  ) noexcept {
    // we lose quite a bit of performance because we are
    // reading from a non changing concurrent map
    // tbb::concurrent_hash_map doesn't support
    // non locking read only access when you can guarantee
    // that the map will not change
    //
    // as observed with graph construction we are ~2x
    // slower on a single thread and we break even with
    // 2 threads - this is acceptable as we don't expect
    // to run this with less than 4 threads, and usually
    // we use at least 32

    de_bruijn_graph_t::const_accessor accessor;
    dbg.find(accessor, starting_minimizer);
    auto const* current_node = &(*accessor);

    // accessor.release() does not do much in terms of performance
    // and even though we know the hash_map is not going to change
    // better to avoid it because surely UB


    simplified_graph_t::mapped_type rv{*current_node};
    
    // = \       / =
    // = = ===== = =
    // = /       \ =
    //   X       X
    while (current_node->second.out_edges.size() == 1) {
      dbg.find(accessor, *current_node->second.out_edges.begin());
      current_node = &(*accessor);
      
      if (current_node->second.fan_in) {
        break;
      }

      rv.push_back(*current_node);
    }

    return rv;
  }

  void unitig_task(
    simplified_graph_t& simplified,
    detail::compact_minimizer const minimizer,
    de_bruijn_graph_t const& dbg,
    ::std::mutex& to_process_mutex,
    ::tbb::task_group& to_process
  ) noexcept {
    if (simplified.count(minimizer)) {
      return;
    }

    auto const& chain = unitig(dbg, minimizer);
    for (auto const& out_edge : chain.back().second.out_edges) {
      // checking if out_edge is in simplified here does not yield any performance
      to_process_mutex.lock();
      to_process.run([&, minimizer = out_edge]{
        unitig_task(simplified, minimizer, dbg, to_process_mutex, to_process);
      });
      to_process_mutex.unlock();
    }

    simplified.insert({minimizer, ::std::move(chain)});
  }

  simplified_graph_t simplify(de_bruijn_graph_t const& dbg) noexcept {
    if (dbg.empty()) {
      return simplified_graph_t{};
    }

    ::std::mutex to_process_mutex;
    ::tbb::task_group to_process;

    simplified_graph_t simplified;

    for (auto const& [minimizer_ref, node] : dbg) {
      if (node.fan_in || !node.last_in.has_value()) {
        to_process.run([&, minimizer = minimizer_ref]{
          unitig_task(simplified, minimizer, dbg, to_process_mutex, to_process);
        });
      }
    }

    to_process.wait();
    return simplified;
  }

  void write_gfa(
    ::std::ostream& out,
    simplified_graph_t const& graph,
    sequence_index_t const& index,
    command_line_options const& opts
  ) noexcept {
    auto get_index = [
      indexes = minimizer_map_t<::std::size_t>{}, index = 0
    ](auto&& key) mutable {
      auto const [iter, is_inserted] = indexes.try_emplace(key, index);
      index += is_inserted;
      return iter.value();
    };

    out << "H\tVN:Z:1.0" << "\n";

    for (auto const& [starting_minimizer, minimizer_node_pairs] : graph) {
      if (!out) {
        break;
      }

      out << "S\t" << get_index(starting_minimizer) << "\t";

      ::std::size_t total_len = 0;

      for (auto const& [current_minimizer, node_data] : minimizer_node_pairs) {

        auto const begin = current_minimizer.minimizer;
        auto const len   = 
          calculate_length(
            current_minimizer.minimizer,
            node_data.out_edges.empty() ? opts.k - 1 : 2,
            opts.l
          ) + (node_data.out_edges.empty() ? opts.l : 0);

        total_len += len;
        
        if (opts.sequences) {
          auto const& read = index(begin->read);
          out.write(read.data() + begin->offset,
            static_cast<::std::streamsize>(len));
        }
      }

      if (!opts.sequences) {
        out << "*";
      }

      out << "\t"
          << "LN:i:" << total_len
          << "\n";
    }

    for (auto const& [starting_minimizer, minimizer_node_pairs] : graph) {
      auto const& [prefix_minimizer, prefix_node] = minimizer_node_pairs.back();

      for (auto const& out_key : prefix_node.out_edges) {
        out << "L\t" << get_index(starting_minimizer) 
            << "\t+\t" << get_index(out_key) << "\t+\t"
            << 0 << "M"
            << "\n";
      }
    }

    out << "# cpp-mdbg de Bruijn minimizer graph"
        << "\n"
        << "# "
        << opts
        << "\n";
    
    if (!out) {
      ::mdbg::terminate("unable to write graph to ", opts.output_prefix);
    }
  }

}
