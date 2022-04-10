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

  void write_gfa(
    ::std::ostream& out,
    de_bruijn_graph_t const& graph,
    sequences_t const& reads,
    bool const write_sequences
  ) noexcept {
    out << "H\tVN:Z:1.0" << "\n";

    // some tools won't allow edges L that use sequences S defined later
    // so we define them all at once
    for (::std::size_t i = 0; i < graph.size(); ++i) {
      auto const [begin, end] = graph[i].first;
      auto const len = (end - 1)->offset - begin->offset;

      out << "S\t" << i << "\t";

      if (write_sequences) {
        auto const& read = *reads[begin->read];
        out.write(
          read.data() + begin->offset, 
          static_cast<::std::streamsize>(len));
      } else {
        out << "*";
      }

      out << "\t"
          << "LN:i:" << len
          << "\n";
    }

    for (::std::size_t i = 0; i < graph.size(); ++i) {
      // calculate the length of N - 1 suffix of prefix
      // in other words the shared part equal to the
      // prefix N - 1 of the node we are connecting to
      auto [prefix_begin, prefix_end] = graph[i].first;

      ++prefix_begin;
      --prefix_end;

      auto const prefix_len = prefix_end->offset - prefix_begin->offset;
      assert(prefix_len > 0);

      for (auto const j : graph[i].second) {
        auto [suffix_begin, suffix_end] = graph[j].first;
        suffix_end -= 2;

        auto const suffix_len = suffix_end->offset - suffix_begin->offset;
        assert(suffix_len > 0);

        // TODO: unsubstantiated
        auto const minimizer_span = ::std::min(prefix_len, suffix_len);

        out << "L\t" << i << "\t+\t" << j << "\t+\t"
            << minimizer_span << "M" 
            << "\n";
      }
    }
  }

}
