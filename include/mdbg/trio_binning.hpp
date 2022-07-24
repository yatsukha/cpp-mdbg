#pragma once

#include <mdbg/io.hpp>
#include <mdbg/opt.hpp>

#include <tsl/robin_map.h>
#include <thread_pool_include.hpp>

#include <cstdint>

namespace mdbg {

  using kmer_counts_t = ::tsl::robin_map<::std::uint64_t, ::std::size_t>;

  kmer_counts_t count_kmers(
    sequences_t const& seqs,
    command_line_options const& opts
  ) noexcept;

  void reduce_to_unique_kmers(
    kmer_counts_t& l,
    kmer_counts_t& r
  ) noexcept;

  // TODO: better namespacing
  ::std::pair<sequences_t, sequences_t> filter_reads(
    sequences_t const& seqs,
    ::std::vector<kmer_counts_t> const& counts,
    command_line_options const& opts
  ) noexcept;

}
