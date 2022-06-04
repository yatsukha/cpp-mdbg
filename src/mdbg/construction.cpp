#include <mdbg/opt.hpp>
#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>

#include <iostream>
#include <cassert>
#include <cstdint>
#include <functional>

namespace mdbg {

  namespace detail {

    inline bool collision(
      compact_minimizer const& l, compact_minimizer const& r,
      ::std::size_t const length
    ) noexcept {
      if (!(l.cached_hash == r.cached_hash)) {
        return false;
      }

      auto l_iter = l.minimizer;
      auto r_iter = r.minimizer;
      
      for (::std::size_t i = 0; i < length; ++i) {
        if ((l_iter++)->minimizer != (r_iter++)->minimizer) {
          return true;
        }
      }

      return false;
    }

  }

  void write_gfa(
    ::std::ostream& out,
    de_bruijn_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept {
    out << "H\tVN:Z:1.0" << "\n";

    minimizer_map_t<::std::size_t> indexes;
    auto const length = opts.k - 1;

    auto get_index = [&indexes, index = 0](auto&& key) mutable -> ::std::size_t {
      auto const [iter, is_inserted] = indexes.try_emplace(key, index);
      index += is_inserted;
      return iter.value();
    };

    for (auto const& [key, unused] : graph) {
      auto const begin = key.minimizer;
      auto const len   = calculate_length(key.minimizer, length, opts.l);
      
      out << "S\t" << get_index(key) << "\t";

      if (opts.sequences) {
        auto const& read = *reads[begin->read];
        out.write(read.data() + begin->offset,
          static_cast<::std::streamsize>(len));
      } else {
        out << "*";
      }

      out << "\t"
          << "LN:i:" << len
          << "\n";
    }

    for (auto const& [key, value] : graph) {
      auto const prefix_len = 
        calculate_length(key.minimizer + 1, length - 1, opts.l);

      for (auto const& out_key : value.out_edges) {
        auto const suffix_len =
          calculate_length(out_key.minimizer, length - 1, opts.l);

        out << "L\t" << get_index(key) 
            << "\t+\t" << get_index(out_key) << "\t+\t"
            << ::std::min(prefix_len, suffix_len) << "M"
            << "\n";
      }
    }

    out << "# cpp-mdbg de Bruijn minimizer graph"
        << "\n"
        << "# "
        << opts
        << "\n";
    
    if (!out) {
      ::mdbg::terminate("unable to write graph to ", opts.output);
    }
  }

  void merge_graphs(
    ::std::vector<de_bruijn_graph_t>::iterator begin,
    ::std::vector<de_bruijn_graph_t>::iterator end
  ) noexcept {
    if (begin == end) {
      return;
    }

    auto& merged = *begin;

    while (++begin != end) {
      auto& other = *begin;
      
      for (auto const& [key, value] : other) {
        auto& data = merged[key];

        data.out_edges.insert(
          data.out_edges.end(), 
          value.out_edges.begin(),
          value.out_edges.end());

        data.read_references.insert(
            data.read_references.end(),
            value.read_references.begin(),
            value.read_references.end());
      }
      
      // force free
      de_bruijn_graph_t{}.swap(other);
    }
  }

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t>::const_iterator begin,
    ::std::vector<read_minimizers_t>::const_iterator end,
    command_line_options const& opts
  ) noexcept {
    auto const overlap_length = opts.k - 1;
    if (opts.k < 2) {
      ::mdbg::terminate("k should be at least 2");
    }

    de_bruijn_graph_t graph;
    ::tsl::robin_map<detail::compact_minimizer, char, detail::compact_minimizer_hash, detail::compact_minimizer_eq> collisions;

    while (begin != end) {
      auto const& read_minimizers = *(begin++);
      if (read_minimizers.size() < opts.k) {
        continue;
      }

      detail::compact_minimizer current_window{
        read_minimizers.begin(),
        {}
      };

      for (::std::size_t i = 0; i < overlap_length; ++i) {
        current_window.cached_hash.advance(read_minimizers[i].minimizer);
      }

      auto current_iter = graph.insert({current_window, {}}).first;
      current_iter.value().read_references.push_back(read_minimizers.begin());
      
      for (::std::size_t i = 1; 
           i < read_minimizers.size() - overlap_length + 1; ++i) {

        ++current_window.minimizer;
        current_window.cached_hash.rotate(
          read_minimizers[i + overlap_length - 1].minimizer,
          read_minimizers[i - 1].minimizer,
          overlap_length
        );

        auto prefix_minimizer = current_iter.key();
        auto& prefix_node = current_iter.value();

        // TODO: consider using just a single read
        prefix_node.out_edges.push_back(current_window);
        prefix_node.fan_out = prefix_node.fan_out ||
          prefix_node.out_edges.front().cached_hash
            != prefix_node.out_edges.back().cached_hash;

        auto const [next_iter, inserted] = graph.insert({current_window, {}});

        (current_iter = next_iter).value().read_references.push_back(
          read_minimizers.begin() + 
            static_cast<detail::minimizer_iter_t::difference_type>(i));

        auto& suffix_node = current_iter.value();
        if (!suffix_node.fan_in && suffix_node.last_in.has_value() 
              && suffix_node.last_in.value().cached_hash 
                  != prefix_minimizer.cached_hash) {
          suffix_node.fan_in = true;
        }

        suffix_node.last_in = prefix_minimizer;
        
        if (opts.check_collisions && !inserted) {
          if (detail::collision(current_iter.key(), current_window, overlap_length)) {
            collisions.insert({current_window, '\0'});
          }
        }
      }
    }

    if (opts.check_collisions) {
      ::std::printf("de Bruijn node collisions: %lu\n", collisions.size());
    }

    return graph;
  }

}
