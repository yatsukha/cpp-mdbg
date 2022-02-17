#include <cassert>
#include <functional>
#include <mdbg/construction.hpp>

namespace mdbg {

  auto& get_end_trie(
    dbg_trie_t& trie,
    minimizer_iter_t begin,
    minimizer_iter_t const end
  ) noexcept {
    auto* curr = &trie;
    while (begin != end) {
      curr = &(*curr)[begin->minimizer];
      ++begin;
    }
    return *curr;
  }

  de_bruijn_graph_t construct(
    ::std::vector<read_minimizers_t> const& minimizers,
    ::std::size_t const k
  ) noexcept {
    assert(k > 1);

    dbg_trie_t trie{};
    de_bruijn_graph_t nodes{};

    auto const get_node = 
      [&trie, &nodes](auto&& begin, auto&& end) {
        auto& trie_node = get_end_trie(trie, begin, end);
        if (!trie_node.opt_element) {
          trie_node.opt_element = nodes.size();
          nodes.push_back({{begin, end}, {}});
        }
        return *trie_node.opt_element;
      };

    auto const delta = static_cast<minimizer_iter_t::difference_type>(k - 1);

    for (auto const& read_minimizers : minimizers) {
      if (read_minimizers.size() < k) {
        continue;
      }

      auto const end = read_minimizers.end() - delta + 1;

      auto iter = read_minimizers.begin();
      auto left = get_node(iter, iter + delta);

      while (++iter != end) {
        auto right = get_node(iter, iter + delta);
        nodes[left].second.push_back(right);
        left = right;
      }
    }

    return nodes;
  }

}