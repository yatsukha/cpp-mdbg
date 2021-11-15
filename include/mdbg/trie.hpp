#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>

namespace mdbg {

  template<typename Element, ::std::size_t NChildren>
  struct trie {
    ::std::array<
      ::std::unique_ptr<trie>, 
      NChildren
    > children;

    auto& operator[](::std::size_t const index) {
      if (!children[index]) {
        children[index] = ::std::make_unique<trie>();
      }
      return *children[index];
    }

    ::std::optional<Element> opt_element;
  };

}
