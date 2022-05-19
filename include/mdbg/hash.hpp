#pragma once

#include <cstdint>
#include <ostream>

namespace mdbg {

  // 128 bit hash value meant for use as a rolling hash
  struct hash128 {
    ::std::uint64_t lower = 0;
    ::std::uint64_t upper = 0;

    friend bool operator==(hash128 const& l, hash128 const& r) noexcept {
      return l.lower == r.lower && l.upper == r.upper;
    }

    friend ::std::ostream& operator<<(::std::ostream& out, hash128 const& h) noexcept {
      return out << h.lower << " " << h.upper;
    }

    void advance(::std::uint64_t const in) noexcept {
      upper <<= 1;
      upper |= lower >> 63;

      lower <<= 1;
      lower ^= in;
    }

    hash128() noexcept = default;

    template<typename Iter>
    hash128(Iter begin, Iter const end) noexcept {
      while (begin != end) {
        advance(*begin++);
      }
    }

    void rotate(
      ::std::uint64_t const in,
      ::std::uint64_t const out,
      ::std::size_t const length
    ) noexcept {
      advance(in);

      if (length < 64) {
        lower ^= out << length;
        // TODO: is this correct?
        upper ^= out >> (64 - length);
      } else if (length < 128) { // shifting by 64 bits or more == UB
        upper ^= out << (length - 64);
      }
    }
  };

}
