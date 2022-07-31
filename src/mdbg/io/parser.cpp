#include <mdbg/io/parser.hpp>
#include <mdbg/io/gzreader.hpp>

namespace mdbg::io {

  void parse_fasta(char const* file, fasta_constumer& consumer) noexcept {
    ::std::string name, sequence;
    gzreader reader{file};

    enum state {
      parsing_name,
      parsing_sequence,
      done,
    } current_state = state::parsing_name;

    auto view = reader.read();

    ++view.first;
    --view.second;

    for (;;) {
      switch (current_state) {
        case state::done:
          return;

        case state::parsing_name: {
          auto eol = 
            reinterpret_cast<char const*>(
              ::std::memchr(view.first, '\n', view.second));

          if (eol != nullptr) {
            name.append(view.first, eol);
            ++eol; // skip '\n'
            view.second -= static_cast<::std::size_t>(eol - view.first);
            view.first = eol;

            current_state = state::parsing_sequence;

            // [[fallthrough]];
          } else {
            name.append(view.first, view.first + view.second);

            if (reader.eof()) {
              ::mdbg::terminate("Unexpected EOF in ", file);
            }
            view = reader.read();

            break;
          }
        }

        case state::parsing_sequence: {
          auto ret = ::std::strpbrk(view.first, "\n>");

          if (ret == nullptr) {
            sequence.append(view.first, view.first + view.second);

            if (reader.eof()) {
              current_state = state::done;
              consumer(::std::exchange(name, {}), ::std::exchange(sequence, {}));
            } else {
              view = reader.read();
            }
          } else {
            sequence.append(view.first, ret);

            if (*ret == '>') {
              current_state = state::parsing_name;
              consumer(::std::exchange(name, {}), ::std::exchange(sequence, {}));
            }

            ++ret;
            view.second -= static_cast<::std::size_t>(ret - view.first);
            view.first = ret;
          }

          break;
        }
      }
    }
  }

}
