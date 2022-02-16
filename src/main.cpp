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
#include <mdbg/construction.hpp>
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>
#include <mdbg/fast_custom_factory.hpp>

::std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0};

int main(int argc, char** argv) {
  if (argc != 2) {
    ::mdbg::terminate("Required arguments: <input reads>");
  }

  if (!::std::filesystem::exists(argv[1])) {
    ::mdbg::terminate("Could not locate given file: ", argv[1]);
  }

  auto timer = ::mdbg::timer{};

  // load the sequences, ty @tbrekalo
  auto parser =
      fast::CreateFastaParser<::std::string>(argv[1]);
  auto const seqs = parser.Parse();

  ::std::cout << "loaded " << seqs.size() << " sequences in "
              << timer.reset_ms() << " ms" << "\n";

  auto constexpr l = 7;
  auto constexpr d = 0.08;

  auto const minimizers = ::mdbg::pick_minimizers(l, d);

  ::std::cout << "picked minimizers in " << timer.reset_ms()
              << " ms" << "\n";

  ::std::vector<::std::vector<::mdbg::detected_minimizer>> detected;
  detected.reserve(seqs.size());

  ::std::size_t sum = 0;
  for (::std::size_t i = 0; i < seqs.size(); ++i) {
    detected.emplace_back(::mdbg::detect_minimizers(*seqs[i], i, minimizers));
    sum += detected.back().size();
  }

  sum /= seqs.size();

  ::std::cout << "average detected minimizers per read: " << sum
              << " in " << timer.reset_ms() << " ms" << "\n";
}
