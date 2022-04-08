#include <mdbg/util.hpp>
#include <mdbg/opt.hpp>

#include <cxxopts.hpp>

#include <filesystem>

namespace mdbg {

  command_line_options command_line_options::parse(int argc, char** argv) noexcept {
    ::cxxopts::Options options("mdbg", "C++ minimizer de Bruijn assembler.");
    options.add_options()
      ("k,kmers", "Length of window of minimizers to use for the de Bruijn graph.",
        ::cxxopts::value<::std::size_t>()->default_value("5"))
      ("l,letters", "Length of the minimizers.",
        ::cxxopts::value<::std::size_t>()->default_value("7"))
      ("d,density", "Density of the universe minimizers.",
        ::cxxopts::value<double>()->default_value("0.008"))
      ("dry-run", "Dry run, do not write to output.",
        ::cxxopts::value<bool>()
          ->default_value("0")
          ->implicit_value("1"))
      ("i,input", "Input reads.", ::cxxopts::value<::std::string>())
      ("o,output", "Output for the graph formatted as GFA.",
        ::cxxopts::value<::std::string>());

    options.parse_positional({"input", "output", ""});
    options.positional_help("input.fastx output.gfa");

    if (argc <= 1) {
      ::mdbg::terminate(options.help());
    }

    command_line_options rv{};

    try {
      auto r = options.parse(argc, argv);

      rv.k = r["k"].as<decltype(rv.k)>();
      rv.l = r["l"].as<decltype(rv.l)>();
      rv.d = r["d"].as<decltype(rv.d)>();

      rv.dry_run = r["dry-run"].as<decltype(rv.dry_run)>();

      rv.input  = r["i"].as<decltype(rv.input)>();
      rv.output = r["o"].as<decltype(rv.output)>();
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
        << ", input=" << ::std::filesystem::absolute(opts.input)
        << ", output=" << ::std::filesystem::absolute(opts.output)
        << ")";
    return out;
  }

}
