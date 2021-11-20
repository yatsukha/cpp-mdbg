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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wparentheses"
#include <fast/fast.hpp>
#include <biosoup/nucleic_acid.hpp>
#pragma GCC diagnostic pop

#include <mdbg/sim/read_sim.hpp>
#include <mdbg/construction.hpp>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

template<typename... Args>
[[noreturn]] void terminate(Args&& ...args) noexcept {
  (..., (::std::cerr << args)) << ::std::endl;
  ::std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
  using namespace ::std::string_literals;

  if (argc != 3) {
    terminate("Required arguments: <input> <gfa output>");
  }

  if (!::std::filesystem::exists(argv[1])) {
    terminate("Could not locate given file: ", argv[1]);
  }

  // load the sequences, ty @tbrekalo
  auto parser =
      fast::CreateFastaParser<::biosoup::NucleicAcid>(argv[1]);
  auto const seqs = parser.Parse();

  ::std::cout << "Loaded sequences." << "\n";
 
  // 10th chromosome expected to be the only sequence
  auto const& ch10 = *seqs[0];
  auto const cutoff = 5'000;

  auto start_time = ::std::chrono::high_resolution_clock::now();

  // simulate reads similar to PacBio HiFi
  auto const reads = ::mdbg::sim::simulate(
    ch10,
    cutoff,
    ::mdbg::sim::configurations::TestConfiguration
  );

  // calculate min and avg coverage
  ::std::vector<::std::int32_t> coverage(cutoff + 1, 0);

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

 
  auto time_taken =
    ::std::chrono::duration_cast<::std::chrono::milliseconds>(
      ::std::chrono::high_resolution_clock::now() - start_time
    );

  
  ::std::cout << "simulated reads in: " << time_taken.count() << " ms "
              << ::std::endl;

  start_time = ::std::chrono::high_resolution_clock::now();

  auto const k = 15;
  auto graph = ::mdbg::construct(reads, k);

  time_taken =
    ::std::chrono::duration_cast<::std::chrono::milliseconds>(
      ::std::chrono::high_resolution_clock::now() - start_time
    );

  auto const print_node = [](auto iter, auto& stream) {
    for (::std::size_t i = 0; i < k - 1; ++i, ++iter) {
      stream << ::biosoup::kNucleotideDecoder[*iter];
    }
  };

  ::std::ofstream out{argv[2], ::std::ios_base::out | ::std::ios_base::trunc};
  if (!out.is_open()) {
    terminate("Could not create/truncate output file ", argv[2]);
  }

  for (::std::size_t i = 0; i < graph.size(); ++i) {
    out << "S\t" << i << "\t";
    print_node(graph[i].first, out);
    out << "\n"; 
  }

  for (::std::size_t i = 0; i < graph.size(); ++i) {
    for (auto const next : graph[i].second) {
      out << "L\t" << i << "\t+\t"
          << next << "\t+\t" 
          << (k - 1) << "M" << "\n";
    } 
  }

  ::std::cout << "dbg graph with k " << k 
              << " and size " << graph.size()
              << "\n";
  ::std::cout << "construction time: " << time_taken.count() << " ms" << "\n";
}
