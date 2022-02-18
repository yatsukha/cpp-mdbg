#pragma once

#include <memory>
#include <string>
#include <vector>

namespace mdbg {

  ::std::vector<::std::unique_ptr<::std::string>> 
  load_sequences(::std::string const&) noexcept;

}
