#include <mdbg/opt.hpp>
#include <mdbg/graph/construction.hpp>
#include <mdbg/util.hpp>

#include <iostream>
#include <cassert>
#include <cstdint>
#include <functional>

namespace mdbg::graph {

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

  void construct(
    de_bruijn_graph_t& graph,
    read_minimizers_t const& read_minimizers,
    command_line_options const& opts
  ) noexcept {
    if (opts.k < 3) {
      ::mdbg::terminate("k should be at least 3");
    }
    auto const overlap_length = opts.k - 1;

    if (read_minimizers.size() < opts.k) {
      return;
    }

    detail::compact_minimizer current_window{
      read_minimizers.begin(),
      {}
    };

    for (::std::size_t i = 0; i < overlap_length; ++i) {
      current_window.cached_hash.advance(read_minimizers[i].minimizer);
    }

    de_bruijn_graph_t::accessor accessor;
    graph.insert(accessor, {current_window, {}});
    
    for (::std::size_t i = 1; 
         i < read_minimizers.size() - overlap_length + 1; ++i) {

      ++current_window.minimizer;
      current_window.cached_hash.rotate(
        read_minimizers[i + overlap_length - 1].minimizer,
        read_minimizers[i - 1].minimizer,
        overlap_length
      );

      auto prefix_minimizer = accessor->first;
      auto& prefix_node = accessor->second;

      prefix_node.out_edges.insert(current_window);

      graph.insert(accessor, {current_window, {}});

      auto& suffix_node = accessor->second;
      if (!suffix_node.fan_in && suffix_node.last_in.has_value() 
            && suffix_node.last_in.value().cached_hash 
                != prefix_minimizer.cached_hash) {
        suffix_node.fan_in = true;
      }

      suffix_node.last_in = prefix_minimizer;
    }
  }

  void construct(
    de_bruijn_graph_t& graph,
    ::std::vector<read_minimizers_t>::const_iterator begin,
    ::std::vector<read_minimizers_t>::const_iterator end,
    command_line_options const& opts
  ) noexcept {
    while (begin != end) {
      construct(graph, *(begin++), opts);
    }
  }

}
