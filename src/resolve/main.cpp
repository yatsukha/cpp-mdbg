#include <cstdio>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>
#include <vector>

#include <tsl/robin_map.h>

#include <mdbg/util.hpp>

namespace rsv {

  struct out_edge {
    ::std::size_t node;
    ::std::size_t overlap_length;
  };

  struct node {
    ::std::vector<::std::size_t> in;
    ::std::vector<out_edge> out;

    ::std::size_t length = 0;
  };

  using graph_t = ::std::vector<node>;

  template<typename K, typename V>
  using unordered_map_t = ::tsl::robin_map<K, V>;

  bool is_prefix(
    ::std::string_view const prefix, 
    ::std::string_view const str,
    ::std::size_t const pos = 0
  ) noexcept {
    if (str.length() - pos < prefix.length()) {
      return false;
    }
    auto iter = str.begin() + pos;
    return ::std::equal(iter, iter + prefix.length(), prefix.begin());
  }

  ::std::vector<::std::string> split(
    ::std::string_view const str,
    ::std::string_view const delim
  ) noexcept {
    ::std::vector<::std::string> rv;

    ::std::size_t start = 0;
    ::std::size_t pos   = 0;

    while (pos + delim.length() <= str.length()) {
      if (is_prefix(delim, str, pos)) {
        if (pos == start) {
          start = ++pos;
          continue;
        }

        rv.emplace_back(str.substr(start, pos - start));
        pos += delim.length();
        start = pos;
      } else {
        ++pos;
      }
    }

    if (start != str.length()) {
      rv.emplace_back(str.substr(start, str.length() - start));
    }

    return rv;
  }

  ::std::size_t calc_size(
    graph_t const& graph,
    ::std::vector<bool>& visited,
    ::std::size_t const current
  ) noexcept {
    visited[current] = true;

    ::std::size_t size = graph[current].length;

    for (auto const [out_node, overlap] : graph[current].out) {
      if (!visited[out_node]) {
        size += calc_size(graph, visited, out_node) - overlap;
      }
    }

    return size;
  }

}

int main(int const argc, char const* const* argv) {
  if (argc != 2) {
    ::mdbg::terminate("Program expects a single argument: <gfa file>.");
  }

  auto const& file = argv[1];
  ::std::ifstream in{file};

  if (!in.is_open()) {
    ::mdbg::terminate("Could not open ", file, " for reading.");
  }

  ::rsv::unordered_map_t<::std::string, ::std::size_t> indices;
  ::std::size_t current_index = 0;

  ::rsv::graph_t graph;

  ::std::size_t line_count = 0;
  ::std::string line;
  while (::std::getline(in, line)) {
    ++line_count;

    if (!line.size()) {
      ::std::fprintf(stderr, "malformed line at %ld", line_count);
      continue;
    }

    if (line.front() == '#') {
      continue;
    }

    auto const tokens = ::rsv::split(line, "\t");

    if (!tokens.size() || !tokens[0].size()) {
      ::std::fprintf(stderr, "ignoring invalid line at %ld\n", line_count);
      continue;
    }

    switch (tokens[0][0]) {
      case 'S':
        {
          indices[tokens[1]] = current_index;
          if (graph.size() <= current_index) {
            graph.resize(current_index + 1);
          }
          graph[current_index].length = ::std::stoull(tokens[3].substr(5));

          ++current_index;
        }
        break;
      case 'L':
        {
          auto const from = tokens[1];
          auto const to   = tokens[3];

          auto const from_iter = indices.find(from);

          if (from_iter == indices.end()) {
            ::mdbg::terminate(
              "Invalid use of segment ", from, " on line ", line_count,
              " without previous definition.");
          }

          auto const to_iter = indices.find(to);

          if (to_iter == indices.end()) {
            ::mdbg::terminate(
              "Invalid use of segment ", to, " on line ", line_count,
              " without previous definition.");
          }

          auto const l = from_iter->second;
          auto const r = to_iter->second;
          auto const max = ::std::max(l, r);

          if (graph.size() <= max) {
            graph.resize(max + 1);
          }

          auto const overlap_len = ::std::stoull(tokens[5]);

          graph[l].out.push_back({r, overlap_len});
          graph[r].in.push_back(l);
        }
        break;
      default:
        ::std::fprintf(stderr, "ignoring line with type %c\n", tokens[0][0]);
    }
  }

  ::std::printf("nodes in graph: %ld\n", graph.size());

  ::std::vector<bool> visited(graph.size(), false);
  ::std::size_t graph_size = 0;

  for (::std::size_t i = 0; i < graph.size(); ++i) {
    if (!visited[i]) {
      graph_size += ::rsv::calc_size(graph, visited, i);
    }

    if (!graph[i].in.size()) {
      ::std::printf("%ld has 0 in edges\n", i);
    }

    if (!graph[i].out.size()) {
      ::std::printf("%ld has 0 out edges\n", i);
    }
  }

  ::std::printf("graph size: %ld\n", graph_size);
}
