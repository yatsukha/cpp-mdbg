#include "mdbg/minimizers.hpp"
#include "mdbg/opt.hpp"
#include <mdbg/io.hpp>
#include <mdbg/io/parser.hpp>
#include <mdbg/util.hpp>

#include <filesystem>
#include <mutex>
#include <tbb/task_group.h>

namespace mdbg {

  ::std::string extension(::std::string const& filename) noexcept {
    auto const base = ::std::filesystem::path{filename}.filename().string();
    auto const period = base.find('.');

    if (period == ::std::string::npos) {
      return "";
    }

    return base.substr(period);
  }

  bool fastq(::std::string const& extension) noexcept {
    return
      extension.find(".fastq") != ::std::string::npos ||
      extension.find(".fq") != ::std::string::npos;
  }

  sequences_t convert(
    ::std::vector<::std::unique_ptr<::std::string>>&& seqs
  ) noexcept {
    sequences_t rv(seqs.size());

    for (::std::size_t i = 0; i < rv.size(); ++i) {
      rv[i] = sequences_t::value_type{seqs[i].release()};
    }

    return rv;
  }

  processed_sequences_t load_sequences(
    ::std::string const& input,
    command_line_options const& opts
  ) noexcept {
    if (!::std::filesystem::exists(input)) {
      ::mdbg::terminate("Could not locate given file: ", input);
    }

    if (fastq(extension(input))) {
      ::mdbg::terminate(
          ".fq and .fq.gz are currently unsupported, "
          "use seqtk to convert to .fa or .fa.gz");
    }

    processed_sequences_t rv;
    ::std::mutex rv_lock;
    ::tbb::task_group tg;

    io::fasta_constumer consumer = 
      [&rv, &tg, &opts, &rv_lock, index = 0](
        auto&&, auto&& seq
      ) mutable {
        rv.first.emplace_back(
          ::std::make_shared<::std::string>(
            ::std::forward<::std::string>(seq)));

        tg.run([&rv, &seq = *rv.first.back(), &opts, read_id = index++, &rv_lock]{
          auto minimizers = ::mdbg::detect_minimizers(seq, read_id, opts);
          rv_lock.lock();
          rv.second.emplace_back(::std::move(minimizers));
          rv_lock.unlock();
        });
      };

  ::mdbg::io::parse_fasta(input.c_str(), consumer);

    tg.wait();

    return rv;
  }

}

