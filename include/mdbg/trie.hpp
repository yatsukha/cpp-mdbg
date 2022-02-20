#pragma once

#include <unordered_map>
#include <cstddef>
#include <memory>
#include <optional>

#include <tsl/robin_map.h>

namespace mdbg {

  template<typename T>
  struct trie {
    using key_type   = T;
    using child_type = ::std::unique_ptr<trie>;

    ::tsl::robin_map<key_type, child_type> children;

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
