#pragma once

#include <memory>
#include <vector>
// TODO: consider using https://github.com/sparsehash/sparsehash
#include <unordered_map>

#include <biosoup_include.hpp>

namespace mdbg {

  struct minimizers {
    // TODO: could be costly
    ::std::unordered_map<::std::string, ::std::uint64_t> to_hash;
    ::std::unordered_map<::std::uint64_t, ::std::string> from_hash;

    ::std::size_t length;
  };

  minimizers pick_minimizers(::std::size_t const l, double d) noexcept;

  struct detected_minimizer {
    ::std::size_t read;
    ::std::size_t offset;

    ::std::uint64_t minimizer;
  };

  ::std::vector<detected_minimizer> detect_minimizers(
    ::std::string const& read,
    ::std::size_t const read_id,
    minimizers const& ms
  ) noexcept;

}
