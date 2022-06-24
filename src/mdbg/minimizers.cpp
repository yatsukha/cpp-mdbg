#include "mdbg/opt.hpp"
#include <mdbg/minimizers.hpp>
#include <mdbg/util.hpp>

#include <nthash_include.hpp>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mdbg {

  ::std::vector<detected_minimizer> detect_minimizers(
    ::std::string const& seq,
    ::std::size_t const read_id,
    command_line_options const& opts
  ) noexcept {
    ::std::vector<detected_minimizer> minimizers;
    ::std::uint64_t hash, rc_hash;

    minimizers.reserve(
      static_cast<::std::size_t>(static_cast<double>(seq.size()) * opts.d));

    ::std::uint64_t const integer_density =
      static_cast<::std::uint64_t>(
        opts.d * 
          static_cast<decltype(opts.d)>(
            ::std::numeric_limits<::std::uint64_t>::max()));

    for (::std::size_t i = 0; i < seq.size() - opts.l + 1; ++i) {
      ::std::uint64_t canonical;

      if (i) {
        canonical = ::NTC64(
          static_cast<unsigned char>(seq[i - 1]), 
          static_cast<unsigned char>(seq[i - 1 + opts.l]),
          static_cast<unsigned>(opts.l), hash, rc_hash);
      } else {
        canonical = ::NTC64(
          seq.data(), static_cast<unsigned>(opts.l),
          hash, rc_hash);
      }

      if (canonical <= integer_density) {
        minimizers.push_back({
          read_id, 
          i,
          canonical
        });
      }
    }

    return minimizers;
  }

}
