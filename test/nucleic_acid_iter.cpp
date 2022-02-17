#include <catch2/catch.hpp>

#include <biosoup_include.hpp>
#include <mdbg/nucleic_acid_iter.hpp>

#include <string>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};
using iter = ::mdbg::nucleic_acid_iter;

TEST_CASE("Base case", "[nucleic_acid_iter]") {
  ::std::string sequence{
    "ATGCTGATCGTCGATGCTGATCT"
    "GTCGTGATGCTGATCGTGTCCCC"
  };
  REQUIRE(sequence.length() > iter::inner_acids);

  ::biosoup::NucleicAcid acid("name", sequence);
  REQUIRE(acid.inflated_len == sequence.length());
  REQUIRE(acid.deflated_data.size() ==
    sequence.length() / iter::inner_acids 
      + ((sequence.length() % iter::inner_acids) != 0));

  iter i{acid, 0};

  for (auto const c : sequence) {
    REQUIRE(c == ::biosoup::kNucleotideDecoder[*i]);
    ++i;
  }

  iter begin{acid, 0};
  iter end{acid, sequence.length()};

  REQUIRE(end - begin == sequence.length());

  auto const c_begin = begin;
  while (begin != end) {
    REQUIRE(::biosoup::kNucleotideDecoder[*begin] == sequence[begin - c_begin]);
    ++begin;
  }

  iter indexed{acid, 0};
  for (::std::size_t i = 0; i < sequence.size(); ++i) {
    REQUIRE(sequence[i] == 
      ::biosoup::kNucleotideDecoder[*(indexed + static_cast<::std::int32_t>(i))]);
  }
}

TEST_CASE("Advanced case", "[nucleic_acid_iter]") {
  ::std::string sequence{
    "TCTTCGTATAAAAACTACACAGAATCATTCTCAACAACTACTTTGTGATGTGTGCGTTCAACTCACAAAGTTTAACCTTT" 
    "CTTTTCATAGAGCAGTTTGGAAACACTCTGTTTGTAAAGCCTGCAAGTGCTTTTTTGGACTTCATTGAGGCCTTCGTTGG" 
    "AAACGGGATTTCTTCATATAATGCTAGACAGAAGAATTCTCAGTCACTTCTTTGTGTTGTGTGTATTCAAGTCACAGAGT" 
    "TGAACCTTCCTTTAGACAGAGCAGTTTTGAAAAATTCTTTCTGTGTAATTTGCAAGTGGAGATTTCAAGCGATTTGAGGC" 
    "TAATCTTTGAAATGGAAATATCTTCGTGTAAAAACTGCACAGAATCATTCTCAGAAACTGCTTTGTCATCTGTGCGTTCA" 
  };

  ::biosoup::NucleicAcid acid("", sequence);

  // we're gonna be here for a while
  for (::std::size_t offset = 0; offset <= sequence.size(); ++offset) {
    for (::std::size_t i = 0; i < sequence.size() - offset + 1; ++i) {
      iter begin{acid, i};
      iter end{acid, i + offset};
      auto const c_b = begin;
      while (begin != end) {
        REQUIRE(::biosoup::kNucleotideDecoder[*begin] == sequence[i + (begin - c_b)]);
        ++begin;
      }
    }
  }
}
