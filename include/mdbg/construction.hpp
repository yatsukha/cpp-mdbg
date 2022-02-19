#include <mdbg/trie.hpp>
#include <mdbg/minimizers.hpp>

#include <functional>
#include <ostream>

namespace mdbg {

  // index pointing to an element in some vector
  using dbg_trie_t = trie<::std::size_t>;

  using minimizer_iter_t = read_minimizers_t::const_iterator;

  // TODO: vector in vector
  // depends on lifetime of detected read minimizers
  using de_bruijn_graph_t = 
    // vector of nodes in a graph, each node is uniquelly identified by an index
    ::std::vector< 
      // a node is a pair
      ::std::pair<
        // of k-minimizers, k = end - begin
        ::std::pair<minimizer_iter_t, minimizer_iter_t>,
        // and it's out edges
        ::std::vector<::std::size_t>>>;

  ::std::ostream& operator<<(::std::ostream&, de_bruijn_graph_t const&) noexcept;

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t> const& minimizers,
    ::std::size_t const k
  ) noexcept;

}
