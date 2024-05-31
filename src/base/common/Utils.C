// $Id$

// for parsing of vectors
#include "boost/version.hpp"
#include "boost/regex.hpp"
#include "boost/tokenizer.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/convenience.hpp"

#include "Utils.h"
#include "Messages.h"
#include "RuntimeException.h"

#include <type_vector.h>
#include <tensor_value.h>

#include <cctype>
#include <iostream>
#include <iomanip>
#include <sstream>
// we can use only ANSI header for Win compatibility
#include <sys/time.h>


using namespace std;


Utils::Timer::Timer(void)
{
  reset();
}


Utils::Timer::~Timer(void)
{
}


void
Utils::Timer::reset(void)
{
  struct timezone tz;
  struct timeval now;
  gettimeofday(&now, &tz);
  _start = now.tv_sec + 1e-6 * now.tv_usec;
}


string
Utils::Timer::elapsed_string(void)
{
  struct timeval now;
  struct timezone tz;
  gettimeofday(&now, &tz);

  double now_s = now.tv_sec + 1e-6 * now.tv_usec;
  double diff = now_s - _start;

  int m = ::floor(diff / 60);
  int h = ::floor(m / 60);
  m %= 60;
  double s = diff - 60 * (m + 60 * h);

  ostringstream os;
  os.precision(2);
  os << h << "h " << m << "m " << s << "s";

  return os.str();
}


Utils::Progress::Progress(const string message, unsigned int max_size,
    unsigned int step_size)
{  
  _message = message;
  _progress_size = max_size;
  _progress_step = step_size;
  _progress_counter = 0;
  //_progress_step = (max_size > 100) ? max_size*2 : 100;
 
  ostringstream os;
  os << _message << " progress ";
  Messages::info(os.str(), false);
  cout << "  0% ..." << flush;
}

Utils::Progress::~Progress(void)
{
  Messages::info(" done.");
}

void
Utils::Progress::progress_message(unsigned int progress_counter)
{
  if ( progress_counter % _progress_step == 0 )
  {
    //ostringstream os;
    //os << "\b\b\b\b\b\b\b\b" << setw(3) <<
    //    static_cast<int>(100 * progress_counter / _progress_size) << "% ...";
    //Messages::info(os.str(), false);
    cout << "\b\b\b\b\b\b\b\b" << setw(3) <<
      static_cast<int>(100 * progress_counter / _progress_size) << "% ..." << flush;
  }
}




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


void
Utils::camel_tokenize(const std::string& input, std::vector<std::string>& tokens)
{

  if (input.size()==0) return;

  std::stringstream ss, tok;
  ss << input;

  char c = ss.get();
  if (!isupper(c)) return;

  tok << c;

  while (ss.good())
  {
    c = ss.get();
    if (c == EOF) break;
    if (isupper(c))
    {
      tokens.push_back(tok.str());
      tok.str("");
    }

    tok<<c;
  }

  tokens.push_back(tok.str());
 
}


string
Utils::dirname(const std::string& file)
{
#if BOOST_VERSION >= 104700
  boost::filesystem::path p(file);
#else
  boost::filesystem::path p(file, boost::filesystem::native);
#endif

  return p.branch_path().string();
}



std::string
Utils::basename(const std::string& file)
{
#if BOOST_VERSION >= 104700
  boost::filesystem::path p(file);
#else
  boost::filesystem::path p(file, boost::filesystem::native);
#endif

  return boost::filesystem::basename(p);
}

std::string
Utils::file_extension(const std::string& file)
{
#if BOOST_VERSION >= 104700
  boost::filesystem::path p(file);
#else
  boost::filesystem::path p(file, boost::filesystem::native);
#endif


  return boost::filesystem::extension(p);
}


std::string
Utils::time_to_string(double seconds)
{
  int h = 0;
  int m = 0;
  double s = seconds;

  if (s >= 60.0)
  {
    m = int(floor(s/60.0));
    s -= m * 60.0;
  }

  if (m >= 60)
  {
    h = int(floor(m/60.0));
    m -= h * 60;
  }

  ostringstream os;
  if (h > 0)
    os << h << " h ";

  if (m > 0)
    os << m << " min ";

  os.setf(ios::fixed);
  os.width(2);
  os.precision(2);
  os << s << " sec";

  return os.str();
}



void
Utils::convert_win32_path_to_posix(std::string& path)
{
  size_t n = path.size();

  if (n > 0)
  {
    for (size_t i = 0; i < n; i++)
      if (path[i] == '\\') path[i] = '/';

    // treat drive name:
    // cygwin mounts drives e.g. "C:\" as
    // "/cygdrive/c"
    size_t pos = path.find(':');
    if (pos != string::npos)
    {
      path.erase(pos, 1);
      path = "/cygdrive/" + path;
    }

    // I do not remember if this had any sense:
    //if (path[0] == '/')
    //  path = '/' + path;
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
    else if ((open != close) &&
             (tmp == open && last_letter != '\\')) opened++;
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



void
Utils::trim(std::string& str, const std::string& chars)
{
  if (chars.size() == 0)
    boost::algorithm::trim(str);
  else
    boost::algorithm::trim_if(str, boost::algorithm::is_any_of(chars));

  /*
  std::string::size_type pos = str.find_last_not_of(chars);
  if(pos != std::string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(chars);
    if(pos != std::string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
  */
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
  skip_whitespace(istr);

  // now we are sure to have at least one component and we have for sure
  // one symbol in the stream

  vec.clear();

  tmp = istr.get();

  // for components with the same start and end symbol, e.g. "
  bool str_start = false;

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

      case '"':
        if (str_start)
        {
          // we found a complete vector component
          boost::algorithm::trim(comp);
          if (!comp.empty())
            vec.push_back(convert<T>(comp));
          comp.clear();
          str_start = false;
        }
        else
        {
          str_start = true;

          comp += get_until_matching_symbol(istr, '"', '"');
          //cerr << comp;
          comp.erase(comp.size() - 1, 1);
          istr.unget();
        }
        //cerr << comp << endl;
        break;

      case ',':
      case ' ':
      case '\t':
      case '\n':
        // get rid of white space and empty elements
        tmp = istr.get();
        while ((tmp != EOF) && (std::isspace(tmp) || (tmp == ',')))
          tmp = istr.get();
        istr.unget();

      case ')':
      case ']':
      case '}':
        // we found a complete vector component
        boost::algorithm::trim(comp);
        if (!comp.empty())
          vec.push_back(convert<T>(comp));
        comp.clear();
        break;

      default:
        comp += tmp;
    }
    tmp = istr.get();
  }

  // if there is only a single word
  boost::algorithm::trim(comp);
  if (!comp.empty())
    vec.push_back(convert<T>(comp));

  // we have to get rid of the closing symbol
  //boost::algorithm::trim(comp);
  //if (delete_closing)
  //{
  //  comp.erase(comp.size() - 1, 1);
  //  boost::algorithm::trim(comp);
  //}
  //vec.push_back(convert<T>(comp));
}


template <>
void
Utils::extract_vector(const string& input, vector<double>& vec)
{
  // We first read it as strings. This is to allow for ranges:
  // (0, 0.1, 0.2:0.1:1, 2-5, 5.5)
  vector<string> vs;
  extract_vector(input, vs);

  if (vs.size() > 0)
  {
    vec.resize(0);
    vec.reserve(vs.size());

    vector<string> tok;
    for (size_t i = 0; i < vs.size(); i++)
    {
      // check for matlab style range a:b or a:s:b
      tokenize(vs[i], tok, ":");
      if (tok.size() == 1)
      {
        // check for range a-b
        tokenize(vs[i], tok, "-");
        if ((tok.size() == 2) &&
            (tok[0].find_first_of("eE") == string::npos))
        {
          double a = convert<double>(tok[0]);
          double b = convert<double>(tok[1]);
          double step = (b > a) ? 1.0 : -1.0;
          for (double x = a; x < b; x += step)
            vec.push_back(x);
          vec.push_back(b);
        }
        else
        {
          vec.push_back(convert<double>(vs[i]));
        }
      }
      else if (tok.size() == 2)
      {
        double a = convert<double>(tok[0]);
        double b = convert<double>(tok[1]);
        double step = (b > a) ? 1 : -1;
        for (double x = a; x < b; x += step)
          vec.push_back(x);
        vec.push_back(b);
      }
      else if (tok.size() == 3)
      {
        double a = convert<double>(tok[0]);
        double step = convert<double>(tok[1]);
        double b = convert<double>(tok[2]);
        
        if (step * (b - a) < 0)
          throw RuntimeException("'" + vs[i] + "' is invalid double range");

        double sign = 1.0;
        if (step < 0) sign = -1.0;

        for (double x = a; sign*x < (sign*b + 1e-3*step); x += step)
          vec.push_back(x);
        vec.push_back(b);
      }
      else
        throw RuntimeException("'" + vs[i] + "' is invalid double range");
    }
  }
}



void
Utils::extract_vector(const string& input, libMesh::TypeVector<double>& vec)
{
  vector<double> v;
  extract_vector<double>(input, v);

  switch (v.size())
  {
    case 1:
      vec(0) = v[0];
      vec(1) = v[0];
      vec(2) = v[0];
      break;

    case 2:
      vec(0) = v[0];
      vec(1) = v[0];
      vec(2) = v[1];
      break;

    case 3:
      vec(0) = v[0];
      vec(1) = v[1];
      vec(2) = v[2];
      break;

    default:
      throw RuntimeException("\'" + input + "\' does not represent a 3-vector.");
      break;
  }
}



void
Utils::extract_tensor(const std::string& input, libMesh::RealTensor& tensor)
{
  string in(input);
  trim(in, "() \t\n\r");

  // first, extract rows
  vector<string> rows;
  tokenize(in, rows, ";\n");

  vector<double> row;

  if (rows.size() == 1) // only diagonal is supplied
  {
    extract_vector<double>(rows[0], row);
    switch (row.size())
    {
      case 1: // isotropic case
        tensor(0,0) = tensor(1,1) = tensor(2,2) = row[0];
        break;

      case 2: // one different principal axis, assumed along z
        tensor(0,0) = tensor(1,1) = row[0];
        tensor(2,2) = row[1];
        break;

      case 6: // symmetric tensor in xx yy zz yz xz xy ordering
        tensor(1,2) = tensor(2,1) = row[3];
        tensor(0,2) = tensor(2,0) = row[4];
        tensor(0,1) = tensor(1,0) = row[5];

      case 3: // complete diagonal given
        tensor(0,0) = row[0];
        tensor(1,1) = row[1];
        tensor(2,2) = row[2];
        break;

      default:
        throw RuntimeException("\'" + input + "\' does not represent a 3-tensor.");
    }
  }
  else if (rows.size() == 3)
  {
    for (size_t i = 0; i < 3; ++i)
    {
      extract_vector<double>(rows[i], row);
      if (row.size() != 3)
        throw RuntimeException("\'" + input + "\' does not represent a 3-tensor.");

      for (size_t j = 0; j < 3; ++j)
        tensor(i, j) = row[j];
    }
  }
  else
    throw RuntimeException("\'" + input + "\' does not represent a 3-tensor.");

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
      // check for range a-b
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
      {
        // check for matlab style range a:b or a:s:b
        tokenize(vs[i], tok, ":");
        if (tok.size() == 1)
        {
          vec.push_back(convert<int>(vs[i]));
        }
        else if (tok.size() == 2)
        {
          int a = convert<int>(tok[0]);
          int b = convert<int>(tok[1]);
          int step = (b > a) ? 1 : -1;
          for (int x = a; x <= b; x += step)
            vec.push_back(x);
        }
        else if (tok.size() == 3)
        {
          int a = convert<int>(tok[0]);
          int step = convert<int>(tok[1]);
          int b = convert<int>(tok[2]);
          if (step * (b - a) < 0)
            throw RuntimeException("'" + vs[i] + "' is invalid integer range");
          for (int x = a; x <= b; x += step)
            vec.push_back(x);
        }
        else
          throw RuntimeException("'" + vs[i] + "' is invalid integer range");
      }
    }
  }
}


//template
//void
//Utils::extract_vector<double>(const string& input, vector<double>& vec);

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
Utils::extract_vector<unsigned short>(const string& input, vector<unsigned short>& vec);

template
void
Utils::extract_vector<char>(const string& input, vector<char>& vec);

template
void
Utils::extract_vector<string>(const string& input, vector<string>& vec);

