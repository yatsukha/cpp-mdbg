#pragma once

#include <mdbg/minimizers.hpp>
#include <mdbg/opt.hpp>

#include <memory>
#include <string>
#include <vector>

namespace mdbg {

  using sequences_t = ::std::vector<::std::shared_ptr<::std::string>>;
  using minimizers_t = ::std::vector<read_minimizers_t>;

  using processed_sequences_t = ::std::pair<sequences_t, minimizers_t>;

  processed_sequences_t load_sequences(
    ::std::string const& input,
    command_line_options const& opts
  ) noexcept;

}
