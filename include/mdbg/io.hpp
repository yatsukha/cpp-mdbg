#pragma once

#include <memory>
#include <string>
#include <vector>

namespace mdbg {

  using sequences_t = ::std::vector<::std::shared_ptr<::std::string>>;

  sequences_t load_sequences(::std::string const&) noexcept;

}
