#pragma once

#include <thread>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <type_traits>
#include <utility>

namespace mdbg {

  template<char const* Format, typename Value>
  struct table_cell;

  template<typename... TableCells>
  struct refreshing_table_display;

  template<char const* Format, typename Value>
  struct refreshing_table_display<table_cell<Format, Value>> {
    static_assert(::std::is_integral_v<Value>);

    ::std::atomic<bool> done = false;
    ::std::atomic<::std::remove_cv_t<Value>> counter = 0;

    template<::std::size_t N>
    void increment() noexcept {
      static_assert(N == 0);
      ++counter;
    }

    bool is_done() const noexcept {
      return done.load();
    }

    void print(FILE* out) const noexcept {
      ::std::fprintf(out, Format, counter.load());
      ::std::fflush(out);
    }
  };

  template<char const* Format, typename Value, typename... Rest>
  struct refreshing_table_display<table_cell<Format, Value>, Rest...>
    : refreshing_table_display<Rest...>
  {
    using parent_type = refreshing_table_display<Rest...>;

    static_assert(::std::is_integral_v<Value>);
    ::std::atomic<::std::remove_cv_t<Value>> counter = 0;

    template<::std::size_t N>
    void increment() noexcept {
      if constexpr (N == 0) {
        ++counter;
      } else {
        parent_type::template increment<N - 1>();
      }
    }

    void print(FILE* out) const noexcept {
      ::std::fprintf(out, Format, counter.load());
      parent_type::print(out);
    }
  };

  template<typename... Cells>
  struct table_printer {
    refreshing_table_display<Cells...> table = {};
    ::std::size_t const sleep_duration;
    FILE* output;
    ::std::thread printer_thread;

    table_printer(table_printer const&) = delete;
    table_printer& operator=(table_printer const&) = delete;

    table_printer(::std::size_t const sleep_duration, FILE* output)
      : sleep_duration(sleep_duration)
      , output(output)
      , printer_thread([this]{
          while (!this->table.is_done()) {
            table.print(this->output);
            ::std::this_thread::sleep_for(
              ::std::chrono::milliseconds{this->sleep_duration});
          }
        })
    {}

    ~table_printer() {
      printer_thread.join();
    }
  };


}
