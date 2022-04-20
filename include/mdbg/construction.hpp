#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>

#include <functional>
#include <ostream>

namespace mdbg {

  namespace detail {

  using minimizer_iter_t = read_minimizers_t::const_iterator;

    struct compact_minimizer {
      minimizer_iter_t minimizer;
      ::std::size_t length;
      ::std::size_t cached_hash;
    };

    struct compact_minimizer_hash {
      using value_type = compact_minimizer;
      ::std::size_t operator()(value_type const& m) const noexcept {
        return m.cached_hash;
      }
    };

    struct compact_minimizer_eq {
      using value_type = compact_minimizer;
      bool operator()(
        value_type const& l,
        value_type const& r
      ) const noexcept {
        // TODO: prime for 128 or 256 simd
        //       but might slowdown rest of the program depending on the cpu arch

        auto const K = l.length;

        if (l.cached_hash != r.cached_hash) {
          return false;
        }

        auto l_iter = l.minimizer;
        auto r_iter = r.minimizer;
        
        // TODO: test this
        for (::std::size_t i = 0; i < K; ++i) {
          if ((l_iter++)->minimizer != (r_iter++)->minimizer) {
            l_iter = l.minimizer;
            r_iter = r.minimizer + static_cast<long>(K) - 1;

            for (::std::size_t j = 0; j < K; ++j) {
              if ((l_iter++)->minimizer != (r_iter--)->minimizer) {
                return false;
              }
            }
            
            break;
          }
        }

        return true;
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

  using de_bruijn_graph_t = 
    ::tsl::robin_map<
      detail::compact_minimizer, 
      detail::dbg_node,
      detail::compact_minimizer_hash,
      detail::compact_minimizer_eq
    >;

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t> const& minimizers,
    ::std::size_t const k
  ) noexcept; 

  void write_gfa(
    ::std::ostream& out,
    de_bruijn_graph_t const& graph,
    sequences_t const& reads,
    command_line_options const& opts
  ) noexcept;

}
