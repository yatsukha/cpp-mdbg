#pragma once

#include <optional>
#include <ostream>
#include <string>
#include <cstdint>

namespace mdbg {

  struct trio_binning_options {
    ::std::string input_0;
    ::std::string input_1;

    ::std::size_t kmer_length;
    ::std::size_t lower_threshold;
  };

  struct command_line_options {
    ::std::size_t threads;
    ::std::size_t k;
    ::std::size_t l;
    double d;

    bool analysis;
    bool dry_run;
    bool sequences;
    // bool check_collisions;

    ::std::optional<trio_binning_options> trio_binning = ::std::nullopt;

    ::std::string input;
    ::std::string output_prefix;

    static command_line_options parse(int argc, char** argv) noexcept;

    friend ::std::ostream& operator<<(
      ::std::ostream&, command_line_options const&) noexcept;
  };

}
