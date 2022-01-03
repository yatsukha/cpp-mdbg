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
#include <mdbg/util.hpp>
#include <mdbg/minimizers.hpp>

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
      fast::CreateFastaParser<::biosoup::NucleicAcid>(argv[1]);
  auto const seqs = parser.Parse();

  ::std::cout << "loaded " << seqs.size() << " sequences in "
              << timer.reset_ms() << " ms" << "\n";

  auto constexpr l = 13;
  auto constexpr d = 0.0008;

  auto const minimizers = ::mdbg::pick_minimizers(l, d);

  ::std::cout << "picked minimizers in " << timer.reset_ms()
              << " ms" << "\n";
}
