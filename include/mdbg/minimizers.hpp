#pragma once

#include <mdbg/opt.hpp>

#include <memory>
#include <vector>

#include <tsl/robin_map.h>
#include <biosoup_include.hpp>

namespace mdbg {

  struct detected_minimizer {
    ::std::size_t read;
    ::std::size_t offset;

    ::std::uint64_t minimizer;
  };

  using read_minimizers_t = ::std::vector<detected_minimizer>;

  read_minimizers_t detect_minimizers(
    ::std::string const& read,
    ::std::size_t const read_id,
    command_line_options const& opts
  ) noexcept;

}
