#include <mdbg/trie.hpp>
#include <mdbg/minimizers.hpp>

#include <functional>

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

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t> const& minimizers,
    ::std::size_t const k
  ) noexcept;

  // ::std::vector<detected_minimizer>
  //using dbg_trie_t = trie<::std::size_t, 1 << 4>;

  //namespace detail {
  //
  //  inline ::std::reference_wrapper<dbg_trie_t> get_end_trie(
  //    dbg_trie_t& trie,
  //    nucleic_acid_iter begin,
  //    nucleic_acid_iter end
  //  ) noexcept {
  //    ::std::reference_wrapper<dbg_trie_t> curr = trie;

  //    while (!(begin == end)) {
  //      auto letter = *begin;

  //      if (++begin == end) {
  //        curr = curr.get()[letter];
  //        break;
  //      }

  //      curr = curr.get()[static_cast<::std::size_t>((letter << 2) | *begin)];
  //      ++begin;
  //    }

  //    return curr;
  //  }

  //}

  //inline auto construct(
  //  // reads to construct from
  //  ::std::vector<sim::read> const& reads, 
  //  // size of k-mers to use for construction
  //  ::std::size_t const k
  //) noexcept {
  //  dbg_trie_t trie;

  //  ::std::vector<::std::pair<
  //    nucleic_acid_iter, ::std::vector<::std::size_t>
  //  >> dbg_nodes;

  //  for (auto const [_, read_begin, read_end] : reads) {
  //    auto begin = read_begin;
  //    auto end   = begin + static_cast<::std::int32_t>(k - 1);

  //    auto current_trie = detail::get_end_trie(trie, begin, end);
  //    if (!current_trie.get().opt_element) {
  //      current_trie.get().opt_element = dbg_nodes.size();
  //      dbg_nodes.push_back({begin, {}});
  //    }

  //    do {
  //      ++begin; ++end;
  //      auto next_trie = detail::get_end_trie(trie, begin, end);

  //      if (!next_trie.get().opt_element) {
  //        next_trie.get().opt_element = dbg_nodes.size();
  //        dbg_nodes.push_back({begin, {}});
  //      }
  //      
  //      dbg_nodes[*current_trie.get().opt_element].second.push_back(
  //        *next_trie.get().opt_element
  //      );

  //      current_trie = next_trie;
  //    } while (!(end == read_end));
  //  }

  //  return dbg_nodes;
  //}

}
