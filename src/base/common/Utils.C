// $Id$

// for parsing of vectors
#include "boost/regex.hpp"
#include "boost/tokenizer.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/convenience.hpp"

#include "Utils.h"


#include <cctype>
#include <iostream>
#include <sstream>


using namespace std;


string
Utils::extract_typename(const type_info& info)
{
  const char* s = info.name();

  if (s[0] == 'P')
    s++;

  while (isdigit(s[0]))
    s++;

  return s;
}




void
Utils::tokenize(const std::string& input, std::vector<std::string>& tokens,
    const char* delimiter)
{
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(delimiter);
  
  tokenizer tok(input, sep);

  tokens.resize(0);

  tokenizer::iterator it = tok.begin();
  const tokenizer::iterator end = tok.end();
  for (unsigned int k = 0 ; it != end; ++it, k++)
    tokens.push_back(*it);
}




string
Utils::dirname(const std::string& file)
{
  boost::filesystem::path p(file, boost::filesystem::native);

  return p.branch_path().string();
}



std::string
Utils::basename(const std::string& file)
{
  boost::filesystem::path p(file, boost::filesystem::native);

  return boost::filesystem::basename(p);
}



void 
Utils::convert_path_to_unix(std::string& path)
{
  size_t n = path.size();

  if (n > 0)
  {
    for (size_t i = 0; i < n; i++)
      if (path[i] == '\\') path[i] = '/';

    if (path[0] == '/')
      path = '/' + path;
  }
}



const std::string
Utils::get_until_matching_symbol(std::istream& istr, char open, char close)
{
  std::string str = ""; 
  int tmp = 0;
  int last_letter = 0;
  int opened = 1;
  while (1)
  {
    last_letter = tmp;
    tmp = istr.get();
    if (tmp == EOF) return str;
    else if (tmp == open && last_letter != '\\') opened++;
    else if (tmp == close && last_letter != '\\')
    {
      opened--;
      // un-backslashed closing symbol => it's the end of the string
      if (opened == 0) return (str + close);
      else if (tmp == '\\' && last_letter != '\\') 
	continue;  // do not append an unbackslashed backslash
    }
    str += tmp;
  }

  return str;
}



void
Utils::skip_whitespace(std::istream& istr)
{
  int tmp = istr.get();
  while (isspace(tmp))
  {
    tmp = istr.get();
    if (!istr) return;
  }
  istr.unget();
}



template <typename T>
void
Utils::extract_vector(const string& input, vector<T>& vec)
{
  istringstream istr(input);
  skip_whitespace(istr);

  // if input is empty we immediately return
  if (istr.str().size() == 0) return;

  char tmp = istr.get();
  char opening = tmp;

  size_t n = input.size() - 1;
  while (n >= 0)
  {
    tmp = input[n];
    if (!isspace(tmp)) break;
    n--;
  }

  // We require the string to begin with an opening brace and
  // to end with the corresponding closing brace
  bool delete_closing = true;
  switch (opening)
  {
    case '(':
      if (tmp != ')') return;
      break;
   
    case '{':
      if (tmp != '}') return;
      break;

    case '[':
      if (tmp != ']') return;
      break;

    default:
      istr.unget();
      delete_closing = false;
  }

  // now we are sure to have at least one component and we have for sure
  // one symbol in the stream

  vec.clear();

  tmp = istr.get();

  string comp;
  while (tmp != EOF)
  {
    switch (tmp)
    {
      case '(':
        comp += '(' + get_until_matching_symbol(istr, '(', ')');
        break;

      case '{':
        comp += '{' + get_until_matching_symbol(istr, '{', '}');
        break;

      case '[':
        comp += '[' + get_until_matching_symbol(istr, '[', ']');
        break;

      case ',':
        // we found a complete vector component
        boost::algorithm::trim(comp);
        vec.push_back(convert<T>(comp));
        comp.clear();
        break;

      default:
        comp += tmp;
    }
    tmp = istr.get();
  }

  // we have to get rid of the closing symbol
  boost::algorithm::trim(comp);
  if (delete_closing)
  {
    comp.erase(comp.size() - 1, 1);
    boost::algorithm::trim(comp);
  }
  vec.push_back(convert<T>(comp));
}


template <>
void
Utils::extract_vector(const string& input, vector<int>& vec)
{
  // We first read it as strings. This is to allow for ranges:
  // (1, 3, 5-7, 9)
  vector<string> vs;
  extract_vector(input, vs);

  if (vs.size() > 0)
  {
    vec.resize(0);
    vec.reserve(vs.size());

    vector<string> tok;
    for (size_t i = 0; i < vs.size(); i++)
    {
      tokenize(vs[i], tok, "-");
      if ((tok.size() == 2) &&
          (tok[0].find_first_of("eE") == string::npos))
      {
        int a = convert<int>(tok[0]);
        int b = convert<int>(tok[1]);
        int step = (b > a) ? 1 : -1;
        for (int x = a; x <= b; x += step)
          vec.push_back(x);
      }
      else
        vec.push_back(convert<int>(vs[i]));
    }
  }
}


template
void
Utils::extract_vector<double>(const string& input, vector<double>& vec);

template
void
Utils::extract_vector<bool>(const string& input, vector<bool>& vec);


template
void
Utils::extract_vector<unsigned int>(const string& input,
    vector<unsigned int>& vec);

template
void
Utils::extract_vector<short>(const string& input, vector<short>& vec);

template
void
Utils::extract_vector<char>(const string& input, vector<char>& vec);

template
void
Utils::extract_vector<string>(const string& input, vector<string>& vec);

