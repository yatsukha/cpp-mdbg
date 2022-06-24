#include <cstdio>
#include <mdbg/simplification.hpp>
#include <mdbg/construction.hpp>

namespace mdbg {

  simplified_graph_t::mapped_type unitig(
    de_bruijn_graph_t const& dbg,
    detail::compact_minimizer const& starting_minimizer
    ) noexcept {
    auto current_node = *dbg.find(starting_minimizer);
    simplified_graph_t::mapped_type rv{current_node};
    
    // = \       / =
    // = = ===== = =
    // = /       \ =
    //   X       X
    while (!current_node.second.fan_out) {
      if (current_node.second.out_edges.empty()) {
        break;
      }

      current_node = *dbg.find(current_node.second.out_edges.front());
      
      if (current_node.second.fan_in) {
        break;
      }

      rv.push_back(current_node);
    }

    return rv;
  }

  simplified_graph_t simplify(de_bruijn_graph_t const& dbg) noexcept {
    if (dbg.empty()) {
      return {};
    }

    ::std::vector<detail::compact_minimizer> to_process;

    for (auto const& [minimizer, node] : dbg) {
      if (node.fan_in || !node.last_in.has_value()) {
        to_process.push_back(minimizer);
      }
    }

    simplified_graph_t simplified;

    while (to_process.size()) {
      auto const minimizer = to_process.back();
      to_process.pop_back();

      if (simplified.count(minimizer)) {
        continue;
      }

      auto const& chain = (simplified[minimizer] = unitig(dbg, minimizer));
      for (auto const& out_edge : chain.back().second.out_edges) {
        to_process.push_back(out_edge);
      }
    }

    return simplified;
  }

  // TODO: code dup
  void write_gfa(
    ::std::ostream& out,
    simplified_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept {
    out << "H\tVN:Z:1.0" << "\n";

    minimizer_map_t<::std::size_t> indexes;

    auto get_index = [&indexes, index = 0](auto&& key) mutable -> ::std::size_t {
      auto const [iter, is_inserted] = indexes.try_emplace(key, index);
      index += is_inserted;
      return iter.value();
    };

    for (auto const& [starting_minimizer, minimizer_node_pairs] : graph) {
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
          auto const& read = *reads[begin->read];
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
