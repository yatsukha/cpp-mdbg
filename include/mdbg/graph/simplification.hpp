#pragma once

#include <mdbg/graph/construction.hpp>

#include <vector>

#include <tsl/robin_map.h>

namespace mdbg::graph {

  using simplified_graph_t =
    concurrent_minimizer_map_t<::std::vector<de_bruijn_graph_t::value_type>>;

  simplified_graph_t simplify(de_bruijn_graph_t const& dbg) noexcept;

  void write_gfa(
    ::std::ostream& out,
    simplified_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept;

}