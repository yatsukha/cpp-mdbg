#pragma once

#include <cstring>
#include <mdbg/util.hpp>

#include <zlib.h>

#include <memory>
#include <cstddef>

namespace mdbg::io {

  // credit: https://github.com/tbrekalo/fast/blob/master/include/fast/gzutil.hpp
  template<::std::size_t Capacity = 1 << 17>
  class gzreader {
    ::std::array<char, Capacity + 1> buffer;
    ::gzFile file_ptr;
    bool done = false;
      
   public:
    using buffer_view = ::std::pair<char const*, ::std::size_t>;

    gzreader(char const* file) noexcept : file_ptr(::gzopen(file, "r")) {
      if (file_ptr == nullptr) {
        ::mdbg::terminate("Unable to gzopen ", file);
      }
      // internal buffer, does not need to eq capacity
      ::gzbuffer(file_ptr, Capacity);
    }

    ~gzreader() noexcept {
      ::gzclose_r(file_ptr);
    }

    buffer_view read() noexcept {
      auto const len = ::gzread(file_ptr, buffer.data(), buffer.size() - 1);
      if (len == -1) {
        ::mdbg::terminate("Failed to gzread.");
      }

      buffer[len] = '\0';
      done = len != buffer.size() - 1;
      return {buffer.data(), len};
    }

    bool eof() const noexcept {
      return done; 
    }
  };

}
