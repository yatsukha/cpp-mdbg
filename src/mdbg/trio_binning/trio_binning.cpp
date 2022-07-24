#include <mdbg/trio_binning/trio_binning.hpp>

#include <biosoup_include.hpp>

#include <cstdint>
#include <mutex>
#include <tbb/task_group.h>

namespace mdbg::trio_binning {

  kmer_counts_t count_kmers(
    sequences_t const& seqs,
    command_line_options const& opts
  ) noexcept {
    kmer_counts_t count;

    auto const k = opts.trio_binning->kmer_length;
    auto const mask = ~(static_cast<::std::uint64_t>(-1) << (2 * k));

    for (auto const& seq : seqs) {
      if (seq->length() < k) {
        continue;
      }

      ::std::uint64_t kmer = 0;
      for (::std::size_t i = 0; i < k - 1; ++i) {
        kmer = (kmer << 2) 
                | ::biosoup::kNucleotideCoder[static_cast<::std::size_t>((*seq)[i])];
      }

      for (::std::size_t i = k - 1; i < seq->length(); ++i) {
        kmer = (kmer << 2) 
                | ::biosoup::kNucleotideCoder[static_cast<::std::size_t>((*seq)[i])];
        kmer &= mask;

        ++count[kmer];
      }
    }

    auto iter = count.begin();

    while (iter != count.end()) {
      if (iter->second < opts.trio_binning->lower_threshold) {
        iter = count.erase(iter);
      } else {
        ++iter;
      }
    }

    return count;
  }

  void reduce_to_unique_kmers(
    kmer_counts_t& l,
    kmer_counts_t& r
  ) noexcept {
    auto& small = l.size() <= r.size() ? l : r;
    auto& large = l.size() <= r.size() ? r : l;

    auto iter = small.begin();

    while (iter != small.end()) {
      if (large.count(iter.key())) {
        large.erase(iter.key());
        iter = small.erase(iter);
      } else {
        ++iter;
      }
    }
  }

  ::std::pair<sequences_t, sequences_t> filter_reads(
    sequences_t const& seqs,
    ::std::vector<kmer_counts_t> const& counts,
    command_line_options const& opts
  ) noexcept {
    // TODO: hardcoded parameters for filtering
    ::std::size_t constexpr min_kmers = 10;
    float constexpr min_ratio = 1.5;
    // TODO: code dup
    
    ::std::pair<sequences_t, sequences_t> rv;
    // TODO: concurrent collection
    ::std::pair<::std::mutex, ::std::mutex> filtered_guards{};

    // TODO: first and second getting annoyting
    auto const add_first = [&rv, &filtered_guards](auto&& seq) {
      ::std::lock_guard<::std::mutex> guard(filtered_guards.first);
      rv.first.push_back(seq);
    };
    auto const add_second = [&rv, &filtered_guards](auto&& seq) {
      ::std::lock_guard<::std::mutex> guard(filtered_guards.second);
      rv.second.push_back(seq);
    };

    auto const k = opts.trio_binning->kmer_length;
    auto const mask = ~(static_cast<::std::uint64_t>(-1) << (2 * k));
    
    ::tbb::task_group task_group;

    for (auto const& seq : seqs) {
      task_group.run([&seq, add_first, add_second, k, mask, &counts]{
        ::std::uint64_t kmer = 0;
        ::std::pair<::std::size_t, ::std::size_t> unique_counts = {0, 0};

        for (::std::size_t i = 0; i < k - 1; ++i) {
          kmer = (kmer << 2) 
                  | ::biosoup::kNucleotideCoder[static_cast<::std::size_t>((*seq)[i])];
        }

        for (::std::size_t i = k - 1; i < seq->length(); ++i) {
          kmer = (kmer << 2) 
                  | ::biosoup::kNucleotideCoder[static_cast<::std::size_t>((*seq)[i])];
          kmer &= mask;

          if (counts[0].count(kmer)) {
            ++unique_counts.first;
          }
          if (counts[1].count(kmer)) {
            ++unique_counts.second;
          }
        }

        ::std::size_t max = unique_counts.first;
        ::std::size_t min = unique_counts.second;
        bool first = true;

        if (unique_counts.first < unique_counts.second) {
          ::std::swap(min, max);
          first = false;
        }

        if (max >= min_kmers 
            && (min == 0
                || static_cast<float>(max) / static_cast<float>(min) >= min_ratio)) {
          if (first) {
            add_first(seq);
          } else {
            add_second(seq);
          }
        } else {
          add_first(seq);
          add_second(seq);
        }
      });
    }

    task_group.wait();
    return rv;
  }

}
