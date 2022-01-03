#pragma once

#include <memory>
#include <vector>
// TODO: consider using https://github.com/sparsehash/sparsehash
#include <unordered_map>

#include <biosoup/nucleic_acid.hpp>

namespace mdbg {

  struct minimizers {
    // TODO: could be costly
    ::std::unordered_map<::std::string, ::std::uint64_t> to_hash;
    ::std::unordered_map<::std::uint64_t, ::std::string> from_hash;
  };

  minimizers pick_minimizers(::std::size_t const l, double d) noexcept;

}
