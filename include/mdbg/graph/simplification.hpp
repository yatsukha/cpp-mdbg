#pragma once

#include <mdbg/graph/construction.hpp>

#include <vector>
#include <functional>

namespace mdbg::graph {

  using simplified_graph_t =
    concurrent_minimizer_map_t<::std::vector<de_bruijn_graph_t::value_type>>;

  simplified_graph_t simplify(de_bruijn_graph_t const& dbg) noexcept;

  using sequence_index_t = ::std::function<::std::string const&(::std::size_t const)>;

  void write_gfa(
    ::std::ostream& out,
    simplified_graph_t const& graph,
    sequence_index_t const& index,
    command_line_options const& opts
  ) noexcept;

}
