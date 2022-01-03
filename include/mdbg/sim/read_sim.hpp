#pragma once

#include <biosoup/nucleic_acid.hpp>

#include <cstdint>
#include <iostream> // TODO:
#include <functional>
#include <random>
#include <vector>

namespace mdbg::sim {
  
  namespace detail {

    using biosoup_sequence_t = decltype(::biosoup::NucleicAcid::deflated_data);
    
    struct nucleic_acid_iter {
      ::std::uint8_t constexpr static inner_acids = 32;

      biosoup_sequence_t::value_type const* data;
      ::std::uint8_t inner_index;

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
          (l.data - r.data) * nucleic_acid_iter::inner_acids +
            static_cast<::std::int32_t>(l.inner_index) -
            static_cast<::std::int32_t>(r.inner_index));
      }

      friend bool operator==(
        nucleic_acid_iter const& l,
        nucleic_acid_iter const& r
      ) noexcept {
        return l.data == r.data && l.inner_index == r.inner_index;
      }
    };

  }

  struct config {
    // expected read length
    ::std::uint32_t read_length;
    // read length deviation
    ::std::uint32_t read_deviation;
    // minimum desired coverage
    // for simplicity, the coverage at the ends of sequence are somewhat ignored  
    ::std::uint16_t min_coverage;
  };

  namespace configurations {
    
    auto inline constexpr PacBioHiFi = config{
      /*.read_length =*/    20'000ul,
      /*.read_deviation =*/ 5'000ul, 
      /*.min_coverage  =*/  30ul,
    };

    auto inline constexpr TestConfiguration = config{
      /*.read_length =*/    200ul,
      /*.read_deviation =*/ 20ul, 
      /*.min_coverage  =*/  30ul,
    };

  }

  struct read {
    // offset form the start of the complete assembly
    ::std::size_t offset;
    // beginning of the simulated read
    detail::nucleic_acid_iter begin;
    detail::nucleic_acid_iter end;
  };

  inline ::std::vector<read> simulate(
    ::biosoup::NucleicAcid const& sequence,
    ::std::size_t const cutoff,
    config const& config
  ) noexcept {
    // size of the given sample to generate reads from
    ::std::size_t const size = cutoff;
    // how much to shift reads left
    // if we didn't do this the coverage would be larger on end than
    // on the beginning of the sample
    ::std::size_t const shift_left = config.read_length / 2;
    // minimal length of the read, reads shorter than this are discarded
    ::std::size_t const min_len = static_cast<::std::size_t>(
      ::std::max<::std::int32_t>(
        static_cast<::std::int32_t>(config.read_deviation),
        static_cast<::std::int32_t>(config.read_length)
          - 2 * static_cast<::std::int32_t>(config.read_deviation)
      )
    );

    ::std::mt19937 mt{::std::random_device{}()};

    ::std::uniform_int_distribution<::std::size_t>
      read_begin_gen(0, size - 1);

    ::std::normal_distribution<float> read_length_gen(
      static_cast<float>(config.read_length), 
      static_cast<float>(config.read_deviation));

    ::std::vector<read> simulated_reads(
      3 * ((size * config.min_coverage) / config.read_length)
    );

    for (auto& read : simulated_reads) {
      for (;;) {
        auto read_begin = read_begin_gen(mt);
        if (read_begin < shift_left) {
          read_begin = 0;
        } else {
          read_begin -= shift_left;
        }
        
        auto read_end = 
          read_begin + static_cast<decltype(read_begin)>(read_length_gen(mt));

        
        if (read_end - read_begin < min_len) {
          continue;
        }

        read_end = ::std::min(read_end, size);

        read = {
          read_begin,
          detail::nucleic_acid_iter{
            sequence.deflated_data.data()
              + static_cast<::std::int32_t>(
                  read_begin / detail::nucleic_acid_iter::inner_acids),
            static_cast<::std::uint8_t>(
              read_begin % detail::nucleic_acid_iter::inner_acids)
          },
          detail::nucleic_acid_iter{
            sequence.deflated_data.data()
              + static_cast<::std::int32_t>(
                  read_end / detail::nucleic_acid_iter::inner_acids),
            static_cast<::std::uint8_t>(
              read_end % detail::nucleic_acid_iter::inner_acids)
          },
        };

        break;
      }
    }

    return simulated_reads;
  }

}
