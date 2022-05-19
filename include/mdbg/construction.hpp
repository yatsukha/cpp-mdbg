#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>
#include <mdbg/hash.hpp>

#include <functional>
#include <ostream>

namespace mdbg {

  namespace detail {

    using minimizer_iter_t = read_minimizers_t::const_iterator;
    
    struct compact_minimizer {
      minimizer_iter_t minimizer;
      ::std::size_t length;
      // TODO: identify with reverse by having two hashes
      //       a regular and reverse one
      //       two k-min-mers are equal if either of the hashes match with other
      ::mdbg::hash128 cached_hash;
    };

    inline bool collision(
        compact_minimizer const& l, compact_minimizer const& r
    ) noexcept {
      if (!(l.cached_hash == r.cached_hash)) {
        return false;
      }

      auto const K = l.length;

      auto l_iter = l.minimizer;
      auto r_iter = r.minimizer;
      
      for (::std::size_t i = 0; i < K; ++i) {
        if ((l_iter++)->minimizer != (r_iter++)->minimizer) {
          return true;
        }
      }

      return false;
    }

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

    inline ::std::size_t gen_id() noexcept {
      ::std::size_t static current_id = 0;
      return current_id++;
    }

    struct dbg_node {
      ::std::vector<compact_minimizer> out_edges;
      // metadata
      ::std::vector<minimizer_iter_t> read_references;
    };

  }

  // TODO: maybe use a dedicated class?
  using de_bruijn_graph_t = 
    ::tsl::robin_map<
      detail::compact_minimizer, 
      detail::dbg_node,
      detail::compact_minimizer_hash,
      detail::compact_minimizer_eq
    >;

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t>::const_iterator begin,
    ::std::vector<read_minimizers_t>::const_iterator end,
    command_line_options const& opts
  ) noexcept;

  void merge_graphs(
    ::std::vector<de_bruijn_graph_t>::iterator begin,
    ::std::vector<de_bruijn_graph_t>::iterator end
  ) noexcept;

  void write_gfa(
    ::std::ostream& out,
    de_bruijn_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept;

}
