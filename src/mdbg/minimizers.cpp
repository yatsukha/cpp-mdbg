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

  void cartesian_product(
    ::std::string& curr,
    ::std::size_t index,
    ::std::function<void(::std::string const&)> f
  ) noexcept {
    if (index == curr.size()) {
      f(curr);
      return;
    }
    char constexpr static arr[] = {'A', 'C', 'T', 'G'};
    for (auto const c : arr) {
      curr[index] = c;
      cartesian_product(curr, index + 1, f);
    }
  }

  char complement(char const c) noexcept {
    switch (c) {
      case 'A': return 'T';
      case 'T': return 'A';
      case 'G': return 'C';
      case 'C': return 'G';
      default: return '\0';
    }
  }

  bool canonical(::std::string const& s) noexcept {
    for (::std::size_t i = 0; i < s.size(); ++i) {
      if (s[i] > complement(s[s.size() - i - 1])) {
        return false;
      }
    }
    return true;
  }

  minimizers pick_minimizers(::std::size_t const l, double d) noexcept {
    auto const total_strings = static_cast<::std::size_t>(::std::pow(4, l) * d);
    
    minimizers m{};

    m.length = l;
    m.from_hash.reserve(total_strings);
    m.to_hash.reserve(total_strings);

    ::std::string buffer(l, '\0');
    cartesian_product(buffer, 0, [&m, d](::std::string const& lmer) {
      // it is faster to check if a lmer is canonical
      // than to calculate both regular and reverse-complement hashes
      if (canonical(lmer)) {
        auto const hash = ::NTF64(lmer.data(), static_cast<unsigned>(lmer.size()));
        auto const density = 
          static_cast<double>(hash) 
            / static_cast<double>(::std::numeric_limits<decltype(hash)>::max());

        if (density <= d) {
          m.to_hash[lmer] = hash;
          m.from_hash[hash] = lmer;
        }
      }
    });

    return m;
  }

  ::std::vector<detected_minimizer> detect_minimizers(
    ::std::string const& seq,
    ::std::size_t const read_id,
    minimizers const& ms
  ) noexcept {
    ::std::vector<detected_minimizer> minimizers;
    ::std::uint64_t hash, rc_hash;

    for (::std::size_t i = 0; i < seq.size() - ms.length + 1; ++i) {
      auto const canonical = i 
        ? ::NTC64(
            static_cast<unsigned char>(seq[i - 1]), 
            static_cast<unsigned char>(seq[i - 1 + ms.length]),
            static_cast<unsigned>(ms.length), hash, rc_hash)
        : ::NTC64(
            seq.data(), static_cast<unsigned>(ms.length),
            hash, rc_hash);

      if (ms.from_hash.count(canonical)) {
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
