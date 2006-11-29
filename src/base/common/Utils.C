// $Id$

#include "Utils.h"

// for parsing of vectors
#include "boost/regex.hpp"
#include "boost/tokenizer.hpp"

#include <cctype>


std::string
Utils::extract_typename(const std::type_info& info)
{
  const char* s = info.name();

  if (s[0] == 'P')
    s++;

  while (isdigit(s[0]))
    s++;

  return s;
}


template <typename T>
void
Utils::extract_vector(const std::string& input, std::vector<T>& vec)
{

  // the regexp to match the vector
  static const boost::regex regexp("[[:space:]]*?(?:(\\()|(\\[)|(\\{)){1}\
(.*)(?(1)\\)|(?(2)\\]|\\})){1}[[:space:]]*?");

  //typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  typedef boost::tokenizer<boost::escaped_list_separator<char> > tokenizer;
  //boost::char_separator<char> sep(",");
  boost::escaped_list_separator<char> sep;

  boost::cmatch matches;

  if (boost::regex_match(input.c_str(), matches, regexp))
  {
    // the regular expression uses 3 subexpressions to identify
    // the type of braces used. Subexpression 0 contains the whole
    // string, so the vector without braces is subexpression 4
    int n = matches.size();
    if (n == 5)
    {
      std::string match(matches[4].first, matches[4].second);

      // cut the matched string into tokens
      tokenizer tokens(match, sep);

      tokenizer::iterator it = tokens.begin();
      const tokenizer::iterator end = tokens.end();

      // we resize the vector and fill in the found values
      vec.resize(0);
      for ( ; it != end; ++it)
        vec.push_back(Utils::convert<T>((*it)));
    }
  }
}


template
void
Utils::extract_vector<double>(const std::string& input, std::vector<double>& vec);

template
void
Utils::extract_vector<bool>(const std::string& input, std::vector<bool>& vec);

template
void
Utils::extract_vector<int>(const std::string& input, std::vector<int>& vec);

template
void
Utils::extract_vector<char>(const std::string& input, std::vector<char>& vec);

template
void
Utils::extract_vector<std::string>(const std::string& input,
    std::vector<std::string>& vec);

