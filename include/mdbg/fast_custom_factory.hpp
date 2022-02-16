#pragma once

#include <fast/seq_factory.hpp>
#include <fast_include.hpp>

#include <memory>
#include <string>

namespace mdbg {

  struct string_fasta_factory {
    
    ::std::unique_ptr<::std::string> CreateSequence(
        char const*, char const*,
        char const* begin, char const* end) const {
      return ::std::make_unique<::std::string>(::std::string(begin, end));
    }

  };

  struct string_fastq_factory {
    auto CreateSequence(
        char const* name_begin, char const* name_end,
        char const* begin, char const* end,
        char const*, char const*) const {
      return string_fasta_factory{}.CreateSequence(name_begin, name_end, begin, end);
    }
  };

}

namespace fast {

  template<>
  struct FastaFactoryFor<::std::string> {
    using Type = ::mdbg::string_fasta_factory; 
    static_assert(IsFastaFactoryForV<Type, ::std::string>);
  };

  template<>
  struct FastqFactoryFor<::std::string> {
    using Type = ::mdbg::string_fastq_factory; 
    static_assert(IsFastqFactoryForV<Type, ::std::string>);
  };

}
