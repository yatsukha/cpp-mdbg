#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>

#include <fast_include.hpp>
#include <biosoup_include.hpp>

#include <mdbg/sim/read_sim.hpp>
#include <mdbg/util.hpp>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

template<typename... Args>
[[noreturn]] void terminate(Args&& ...args) noexcept {
  (..., (::std::cerr << args)) << ::std::endl;
  ::std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
  using namespace ::std::string_literals;

  if (argc != 4) {
    terminate("Required arguments: <input> <simulated reads output> <reference output>");
  }

  if (!::std::filesystem::exists(argv[1])) {
    terminate("Could not locate given file: ", argv[1]);
  }

  auto timer = ::mdbg::timer{};

  // load the sequences, ty @tbrekalo
  auto parser =
      fast::CreateFastaParser<::biosoup::NucleicAcid>(argv[1]);
  auto const seqs = parser.Parse();

  ::std::cout << "loaded sequences in: "
              << timer.reset_ms() << " ms" << "\n";

  // 10th chromosome expected to be the only sequence
  auto const start_cutoff = 50'000'000u;
  auto const section_length = 50'000'000u;

  auto const ch10 = 
    ::biosoup::NucleicAcid{"", seqs[0]->InflateData().substr(start_cutoff, section_length)};

  ::std::cout << "size of ch10: " << ch10.inflated_len << "\n";
  ::std::cout << "dumping reference" << "\n";

  ::std::ofstream{argv[3]} 
    << ">ch10 start" << start_cutoff << " len " << section_length << "\n" 
    << ch10.InflateData();

  // simulate reads similar to PacBio HiFi
  auto const reads = ::mdbg::sim::simulate(
    ch10,
    section_length,
    ::mdbg::sim::configurations::PacBioHiFi
  );

  // calculate min and avg coverage
  ::std::vector<::std::int32_t> coverage(section_length + 1, 0);

  for (auto const& r : reads) {
    ++coverage[r.offset];
    --coverage[r.offset + (r.end - r.begin)];
  }

  ::std::size_t min_cov = ::std::numeric_limits<decltype(min_cov)>::max();
  ::std::size_t total_cov = 0;
  ::std::size_t curr_cov = 0;

  for (::std::size_t c = 0; c < coverage.size() - 1; ++c) {
    curr_cov += static_cast<decltype(curr_cov)>(coverage[c]);
    min_cov = ::std::min(min_cov, curr_cov);
    total_cov += curr_cov;
  }

  ::std::cout << "min coverage: " << min_cov << "\n"
              << "avg coverage: " << (total_cov / (coverage.size() - 1)) << "\n";
  
  ::std::cout << "simulated reads in: " << timer.reset_ms() << " ms "
              << ::std::endl;

  // output to sim file
  ::std::ofstream out{argv[2]};
  if (!out.is_open()) {
    terminate("Could not open/create output file: ", argv[2]);
  }

  for (auto [offset, begin, end] : reads) {
    out << ">" << ch10.name << "|" << offset << "|" << ::std::to_string(end - begin);
    out << "\n";
    while (begin != end) {
      out << ::biosoup::kNucleotideDecoder[*begin];
      ++begin;
    }
    out << "\n";
  }

  ::std::cout << "written simulated reads to " << argv[2] 
              << " in  " << timer.reset_ms() << " ms " << "\n";
}
