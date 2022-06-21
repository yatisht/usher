#ifndef TEXT_PARSER_H
#define TEXT_PARSER_H

#include <string>
#include <string_view>

class text_parser {
  public:
    text_parser(const std::string &path);
    ~text_parser();

    std::string_view line() const { return std::string_view(curr, eol - curr); }

    size_t get_size() const { return size; }

    void next_line();

    bool done() const { return curr == end; }


    std::string_view get_value(int column) const {
        auto offset = curr;
        auto num_tabs = 0;
        for (auto p = curr; p != eol; ++p) {
            if (*p == '\t') {
                if (column == num_tabs) {
                    return std::string_view(offset, p - offset);
                }
                ++num_tabs;
                offset = p + 1; // Works for chars only
            }
        }
        return std::string_view(offset, eol - offset);
    }


  private:
    size_t size;
    const char *data, *curr, *end, *eol;
};

#endif
