#ifndef _FTK_STRING_HH
#define _FTK_STRING_HH

#include <ftk/ftk_config.hh>
#include <string>

namespace ftk {

static inline bool ends_with(std::string const & value, std::string const & ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// https://stackoverflow.com/questions/9435385/split-a-string-using-c11
static inline std::vector<std::string> split(const std::string& input, const std::string& regex) {
    // passing -1 as the submatch index parameter performs splitting
    std::regex re(regex);
    std::sregex_token_iterator
        first{input.begin(), input.end(), re, -1},
        last;
    return {first, last};
}

}

#endif
