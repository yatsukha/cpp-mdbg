#include <mdbg/util.hpp>
#include <mdbg/opt.hpp>

#include <cxxopts.hpp>

#include <filesystem>
#include <stdexcept>

namespace mdbg {

  // TODO: simplify a bit
  trio_binning_options parse_trio_binning(::std::string const& arg) noexcept {
    auto const assert_not_npos = [](::std::size_t const index) {
      if (index == ::std::string::npos) {
        ::mdbg::terminate("Expected ':', found end of argument for trio binning.");
      }
    };

    trio_binning_options rv{};

    auto separator_index = arg.find(":");
    assert_not_npos(separator_index);

    try {
      rv.kmer_length = ::std::stoul(arg.substr(0, separator_index));

      if (rv.kmer_length > 32) {
        ::mdbg::terminate("Trio binning K must be <= 32.");
      }

      auto const start = separator_index + 1;

      separator_index = arg.find(":", start);
      assert_not_npos(separator_index);

      rv.lower_threshold = ::std::stoul(arg.substr(start, separator_index - start));
    } catch (::std::logic_error const&) {
      ::mdbg::terminate("Expected non-negative number for K/T in trio binning argument.");
    }

    auto const start = separator_index + 1;
    separator_index = arg.find(":", start);
    assert_not_npos(separator_index);

    rv.input_0 = arg.substr(start, separator_index - start);
    rv.input_1 = arg.substr(separator_index + 1);

    return rv;
  }

  command_line_options command_line_options::parse(int argc, char** argv) noexcept {
    ::cxxopts::Options options(
        "mdbg", 
        "C++ minimizer based de Bruijn assembler.");

    options.add_options()
      ("t,threads", 
        "Maximum number of concurrent threads. "
        "NOTE: Default of 0 means max concurrency.",
        ::cxxopts::value<::std::size_t>()->default_value("0"))
      ("k,kmers", "Length of window of minimizers to use for the de Bruijn graph.",
        ::cxxopts::value<::std::size_t>()->default_value("5"))
      ("l,letters", "Length of the minimizers.",
        ::cxxopts::value<::std::size_t>()->default_value("7"))
      ("d,density", "Density of the universe minimizers.",
        ::cxxopts::value<double>()->default_value("0.008"))
      ("a,analysis",
        "Exit after outputting minimizer statistics for given reads.",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("dry-run", "Dry run, do not write to output.",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("s,sequences", 
        "Output sequences contained within minimizers in output GFA.",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("u,unitigs",
       "Simplify straight portions of the graph into unitigs. "
       "Incurs additional running time!",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("c,check-collisions",
       "Check for node collisions when building the de Bruijn graph. "
       "Incurs runtime overhead!",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("i,input", "Input reads.", ::cxxopts::value<::std::string>())
      ("o,output", "Output prefix for the graph(s) formatted as GFA.",
        ::cxxopts::value<::std::string>())
      ("trio-binning",
        "Format: K:T:reads0.fa:reads1.fa\n"
        "Enables trio binning using K length kmers for counting; "
        "discards kmers with frequency below T.\n"
        "K must be <= 32.",
        ::cxxopts::value<::std::string>()
          ->default_value(""));

    options.parse_positional({"input", "output", ""});
    options.positional_help("input.fastx output.gfa");

    if (argc <= 1) {
      ::mdbg::terminate(options.help());
    }

    command_line_options rv{};

    try {
      auto r = options.parse(argc, argv);

      rv.threads = r["threads"].as<decltype(rv.threads)>();
      rv.k = r["k"].as<decltype(rv.k)>();
      rv.l = r["l"].as<decltype(rv.l)>();
      rv.d = r["d"].as<decltype(rv.d)>();

      rv.analysis = r["analysis"].as<decltype(rv.analysis)>();
      rv.dry_run = r["dry-run"].as<decltype(rv.dry_run)>();
      rv.unitigs = r["unitigs"].as<decltype(rv.unitigs)>();
      rv.sequences = r["sequences"].as<decltype(rv.sequences)>();
      rv.check_collisions = r["check-collisions"].as<decltype(rv.check_collisions)>();

      if (auto const& trio_binning_arg = r["trio-binning"].as<::std::string>();
          trio_binning_arg.length() > 0) {
        rv.trio_binning = parse_trio_binning(trio_binning_arg);
      }

      rv.input  = r["i"].as<decltype(rv.input)>();
      rv.output_prefix = r["o"].as<decltype(rv.output_prefix)>();
    } catch (::cxxopts::OptionException const& exc) {
      ::mdbg::terminate(exc.what());
    }

    return rv;
  }

  ::std::ostream& operator<<(
    ::std::ostream& out, command_line_options const& opts
  ) noexcept {
    out << "command_line_options(k=" << opts.k
        << ", l=" << opts.l
        << ", d=" << opts.d
        << ", unitigs=" << opts.unitigs
        << ", sequences=" << opts.sequences
        << ", input=" << ::std::filesystem::absolute(opts.input)
        << ", output=" << ::std::filesystem::absolute(opts.output_prefix)
        << ", trio-binning=";

    if (opts.trio_binning.has_value()) {
      auto const& tb = opts.trio_binning.value();
      out << "(k=" << tb.kmer_length
          << ", t=" << tb.lower_threshold
          << ", reads0=" << tb.input_0
          << ", reads1=" << tb.input_1
          << ")";
    } else {
      out << "none";
    }

    return out << ")";
  }

}
