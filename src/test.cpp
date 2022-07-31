#include <iostream>
#include <mdbg/io/parser.hpp>

int main(int const, char const* const* argv) {
  ::mdbg::io::fasta_constumer lambda = [](auto&& name, auto&& seq) {
    ::std::cout << ">" << name << "\n" << seq << "\n";
  };
  ::mdbg::io::parse_fasta(argv[1], lambda);
}
