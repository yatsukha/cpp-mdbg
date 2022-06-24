#pragma once

#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>
#include <mdbg/hash.hpp>

#include <functional>
#include <optional>
#include <ostream>

namespace mdbg {

  // TODO: rethink having this in detail
  namespace detail {

    using minimizer_iter_t = read_minimizers_t::const_iterator;
    
    struct compact_minimizer {
      minimizer_iter_t minimizer;
      // TODO: identify with reverse by having two hashes
      //       a regular and reverse one
      //       two k-min-mers are equal if either of the hashes match with other
      ::mdbg::hash128 cached_hash;
    };

    struct compact_minimizer_hash {
      using value_type = compact_minimizer;
      ::std::size_t operator()(value_type const& m) const noexcept {
        return m.cached_hash.collapse();
      }
    };

    struct compact_minimizer_eq {
      using value_type = compact_minimizer;
      bool operator()(
        value_type const& l,
        value_type const& r
      ) const noexcept {
        return l.cached_hash == r.cached_hash;
      }
    };

    struct dbg_node {
      ::std::vector<compact_minimizer> out_edges;
      // metadata
      ::std::vector<minimizer_iter_t> read_references;

      ::std::optional<compact_minimizer> last_in = ::std::nullopt;
      bool fan_in = false;
      bool fan_out = false;
    };

  }

  template<typename V>
  using minimizer_map_t =
    ::tsl::robin_map<
      detail::compact_minimizer,
      V,
      detail::compact_minimizer_hash,
      detail::compact_minimizer_eq
    >;

  using de_bruijn_graph_t = minimizer_map_t<detail::dbg_node>;

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t>::const_iterator begin,
    ::std::vector<read_minimizers_t>::const_iterator end,
    command_line_options const& opts
  ) noexcept;

  void merge_graphs(
    ::std::vector<de_bruijn_graph_t>::iterator begin,
    ::std::vector<de_bruijn_graph_t>::iterator end
  ) noexcept;

  inline ::std::size_t calculate_length(
    decltype(detail::compact_minimizer::minimizer) begin,
    ::std::size_t const length,
    ::std::size_t const
  ) noexcept {
      // calculate the length in bases for
      // 'length' minimizers, without including the last minimizer's
      // own length l
      auto const end = begin + static_cast<::std::int64_t>(length) - 1;
      return end->offset - begin->offset;
  }

  void write_gfa(
    ::std::ostream& out,
    de_bruijn_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept;

}
