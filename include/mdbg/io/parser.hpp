#pragma once

#include <string>
#include <functional>

namespace mdbg::io {

  using fasta_constumer = ::std::function<void(::std::string&&, ::std::string&&)>;

  void parse_fasta(char const* file, fasta_constumer& consumer) noexcept;

}
