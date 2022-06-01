#pragma once

#include <ostream>
#include <string>
#include <cstdint>

namespace mdbg {

  struct command_line_options {
    ::std::size_t threads;
    ::std::size_t k;
    ::std::size_t l;
    double d;

    bool analysis;
    bool dry_run;
    bool unitigs;
    bool sequences;
    bool check_collisions;

    ::std::string input;
    ::std::string output;

    static command_line_options parse(int argc, char** argv) noexcept;

    friend ::std::ostream& operator<<(
      ::std::ostream&, command_line_options const&) noexcept;
  };

}
