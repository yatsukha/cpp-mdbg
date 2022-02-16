#pragma once

#include <mdbg/nucleic_acid_iter.hpp>

#include <biosoup_include.hpp>

#include <cstdint>
#include <iostream> // TODO:
#include <functional>
#include <random>
#include <vector>

namespace mdbg::sim {

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
    nucleic_acid_iter begin;
    nucleic_acid_iter end;
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
          read_begin 
            + static_cast<decltype(read_begin)>(
              ::std::max(read_length_gen(mt), 0.f));

        
        if (read_end - read_begin < min_len) {
          continue;
        }

        read_end = ::std::min(read_end, size);

        read = {
          read_begin,
          nucleic_acid_iter{sequence, read_begin}, 
          nucleic_acid_iter{sequence, read_end}
        };

        break;
      }
    }

    return simulated_reads;
  }

}
