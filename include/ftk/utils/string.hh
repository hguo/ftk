#ifndef _FTK_STRING_HH
#define _FTK_STRING_HH

#include <ftk/config.hh>
#include <string>
#include <algorithm>
#include <regex>
#include <cctype>

namespace ftk {

static inline bool starts_with(std::string const & value, std::string const & starting)
{
  if (value.find(starting) == 0) return true;
  else return false;
}

static inline bool ends_with(std::string const & value, std::string const & ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static inline std::string to_lower_cases(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) {
      return std::tolower(c);
  });
  return str;
}

static inline bool ends_with_lower(const std::string& str, const std::string& ending)
{
  return ends_with(to_lower_cases(str), to_lower_cases(ending));
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
