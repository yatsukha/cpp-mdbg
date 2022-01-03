#include <mdbg/minimizers.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wparentheses"
#include <ntHash/nthash.hpp>
#pragma GCC diagnostic pop

#include <iostream>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mdbg {

  void cartesian_product(
    // TODO: consider using something for predictable allocations
    ::std::vector<::std::string>& dst,
    ::std::string& curr,
    ::std::size_t index,
    ::std::size_t limit
  ) noexcept {
    if (index == limit) {
      dst.push_back(curr);
      return;
    }
    char constexpr static arr[] = {'A', 'C', 'T', 'G'};
    for (auto const c : arr) {
      curr[index] = c;
      cartesian_product(dst, curr, index + 1, limit);
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
    ::std::vector<::std::string> lmers;
    auto const total_strings = static_cast<::std::size_t>(::std::pow(4, l));
    lmers.reserve(total_strings);
    ::std::string buffer(l, '\0');

    if (total_strings * l >= (1 << 30)) {
      ::std::cout << "required GBs for minimizers >= "
                  << static_cast<double>(total_strings * l) / (1 << 30)
                  << "\n";
    }

    cartesian_product(lmers, buffer, 0, l);
    
    minimizers m{};
    for (auto const& lmer : lmers) {
      if (canonical(lmer)) {
        auto const hash = ::NTF64(lmer.data(), static_cast<unsigned>(lmer.size()));
        auto const density = 
          static_cast<double>(hash) 
            / static_cast<double>(::std::numeric_limits<::std::uint64_t>::max());

        if (density <= d) {
          m.to_hash[lmer] = hash;
          m.from_hash[hash] = lmer;
        }
      }
    }

    return m;
  }

}
