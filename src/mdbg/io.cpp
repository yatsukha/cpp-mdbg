#include <mdbg/io.hpp>
#include <mdbg/util.hpp>

#include <fast_include.hpp>

#include <filesystem>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

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

namespace mdbg {

  ::std::vector<::std::unique_ptr<::std::string>> 
  load_sequences(::std::string const& input) noexcept {
    if (!::std::filesystem::exists(input)) {
      ::mdbg::terminate("Could not locate given file: ", input);
    }
    return ::fast::CreateFastaParser<::std::string>(input.c_str()).Parse();
  }

}

