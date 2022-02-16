#pragma once

#include <iostream>
#include <memory>
#include <chrono>
#include <functional>

namespace mdbg {

  template<typename... Args>
  [[noreturn]] void terminate(Args&& ...args) noexcept {
    (..., (::std::cerr << args)) << ::std::endl;
    ::std::exit(EXIT_FAILURE);
  }
  
  namespace detail {
    template<typename T>
    using unique_deleter = ::std::unique_ptr<T, ::std::function<void(T*)>>;
  }

  using scope_guard = detail::unique_deleter<int>;

  template<typename F>
  [[nodiscard]] scope_guard mk_scope_guard(F&& f) noexcept {
    static int i = 0;
    return scope_guard{&i, [f](auto*) { f(); }};
  }

  template<typename F>
  [[nodiscard]] auto defer(F&& f) noexcept {
    return mk_scope_guard(::std::forward<F>(f));
  }

  class timer {
    using time_point = ::std::chrono::time_point<::std::chrono::system_clock>;

    static time_point get_now() noexcept {
      return ::std::chrono::high_resolution_clock::now();
    }

    time_point now = get_now(); 

   public:
    timer() noexcept = default;

    auto get() const noexcept {
      return get_now() - now;
    }

    auto get_ms() const noexcept {
      return ::std::chrono::duration_cast<::std::chrono::milliseconds>(get()).count();
    }

    auto reset() noexcept {
      auto rv = now;
      now = get_now();
      return now - rv;
    }

    auto reset_ms() noexcept {
      return ::std::chrono::duration_cast<::std::chrono::milliseconds>(reset()).count();
    }
  };

}
