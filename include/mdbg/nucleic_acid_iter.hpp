#pragma once

#include "biosoup/nucleic_acid.hpp"
#include <biosoup_include.hpp>

#include <cstdint>
#include <functional>
#include <random>
#include <vector>

namespace mdbg {
  using biosoup_sequence_t = decltype(::biosoup::NucleicAcid::deflated_data);
  
  struct nucleic_acid_iter {
    ::std::uint8_t constexpr static inner_acids = 32;

    biosoup_sequence_t::value_type const* data = nullptr;
    ::std::uint8_t inner_index = 0;

    nucleic_acid_iter() noexcept = default;

    nucleic_acid_iter(
      ::biosoup::NucleicAcid const& acid, ::std::size_t offset
    ) noexcept
      : data(acid.deflated_data.data() + offset / inner_acids)
      , inner_index(offset % inner_acids)
    {}

    nucleic_acid_iter(nucleic_acid_iter const&) noexcept = default;
    nucleic_acid_iter& operator=(nucleic_acid_iter const&) noexcept = default;

    ::std::uint8_t operator*() noexcept {
      return ((*data) >> (inner_index * 2)) & 0b11;
    }

    nucleic_acid_iter& operator++() noexcept {
      if (++inner_index >= inner_acids) /*[[unlikely]]*/ {
        ++data;
        inner_index = 0;
      }
      return *this;
    }

    nucleic_acid_iter& operator+=(::std::int32_t diff) noexcept {
      auto const ptr_skips = diff / inner_acids;
      data += ptr_skips;
      diff -= inner_acids * ptr_skips;

      if (diff > 0) {
        inner_index += static_cast<decltype(inner_index)>(diff);
        if (inner_index >= inner_acids) {
          ++data;
          inner_index -= inner_acids;
        }
      } else {
        auto new_inner = static_cast<::std::int32_t>(inner_index) + diff;
        if (new_inner < 0) {
          --data;
          new_inner += inner_acids; 
        }
        inner_index = static_cast<decltype(inner_index)>(new_inner);
      }
      return *this;
    }

    friend nucleic_acid_iter operator+(
      nucleic_acid_iter const& self, ::std::int32_t diff
    ) noexcept {
      nucleic_acid_iter copy = self;
      copy += diff;
      return copy;
    }

    friend ::std::size_t operator-(
      nucleic_acid_iter const& l,
      nucleic_acid_iter const& r
    ) noexcept {
      return static_cast<::std::size_t>(
        static_cast<::std::int32_t>(l.data - r.data) 
          * nucleic_acid_iter::inner_acids 
          + static_cast<::std::int32_t>(l.inner_index) 
            - static_cast<::std::int32_t>(r.inner_index));
    }

    friend bool operator==(
      nucleic_acid_iter const& l,
      nucleic_acid_iter const& r
    ) noexcept {
      return l.data == r.data && l.inner_index == r.inner_index;
    }

    friend bool operator!=(
      nucleic_acid_iter const& l,
      nucleic_acid_iter const& r
    ) noexcept {
      return !(l == r);
    }
  };

}
