#pragma once

#include <unordered_map>
#include <cstddef>
#include <memory>
#include <optional>

namespace mdbg {

  template<typename T>
  struct trie {
    // TODO: unordered_map is a really bad hash map implementation
    // either add better allocator, or different implementation:
    //  - https://martin.ankerl.com/2019/04/01/hashmap-benchmarks-01-overview/
    using key_type   = T;
    using child_type = ::std::unique_ptr<trie>;

    ::std::unordered_map<key_type, child_type> children;

    auto& operator[](key_type const& index) {
      if (!children.count(index)) {
        children[index] = ::std::make_unique<trie>();
      }
      return *children[index];
    }
    
    // TODO: ergonomy
    ::std::optional<key_type> opt_element;
  };

}
