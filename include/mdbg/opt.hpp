#pragma once

#include <string>
#include <cstdint>

namespace mdbg {

  struct command_line_options {
    ::std::size_t k;
    ::std::size_t l;
    double d;

    ::std::string input;
    ::std::string output;

    static command_line_options parse(int argc, char** argv) noexcept;
  };

}
