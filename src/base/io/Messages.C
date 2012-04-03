// $Id$

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include "Messages.h"
#include "TeeStream.h"
#include "TypeDefs.h"
#include "InitFailedException.h"
#include "TiberCad.h"
#include "Utils.h"
#include "tiber_config.h"

#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace
{
#ifndef _WIN32
  const char* red     = "\033[1;31m";
  const char* yellow  = "\033[1;33m";
  const char* redb    = "\033[1;41m";
  const char* yellowb = "\033[1;43m";
  const char* blue    = "\033[1;34m";
  const char* white   = "\033[1;37m";
  const char* normal  = "\033[0m";
#else
  const char* red     = "";
  const char* yellow  = "";
  const char* redb    = "";
  const char* yellowb = "";
  const char* blue    = "";
  const char* white   = "";
  const char* normal  = "";
#endif
}

using namespace std;


const string
Messages::endl = "\n";

const int
Messages::_max_width = 80;

int
Messages::_indent = 0;

const int
Messages::_indent_width = 2;

int
Messages::_error_count = 0;

int
Messages::_warning_count = 0;

bool
Messages::_interactive = false;

ofstream
Messages::_log;



void
Messages::print_statistics(void)
{
  ostringstream os;
  os << _warning_count << " warnings during run\n"
      << _error_count << " errors detected";
  newline();
  info(os.str());
  newline();
}




void
Messages::set_log_file(const string& logfile)
{
  using namespace boost::filesystem;

  path logpath(logfile);//, native);
#if defined(_WIN32)
  logpath = logpath.branch_path();
#else
  logpath.remove_leaf();
  //logpath.remove_filename();
#endif
  if (logpath.string().size() > 0)
  {
    if (!exists(logpath))
    {
      // we catch any error here without doing anything yet
      try {
        create_directories(logpath);
      }
      catch (...) {}
    }

    if (!(exists(logpath) && is_directory(logpath)))
    {
      string msg("Cannot create or use directory '");
      msg += logpath.string() + "' for logging.";
      throw InitFailedException(msg);
    }
  }

  _log.open(logfile.c_str());
  if (_log.fail() || !_log.good())
    throw InitFailedException("cannot open logfile for writing.");

  time_t now;
  time(&now);
  _log << Messages::endl;
  _log << "Log start: " << ctime(&now) << Messages::endl;
  _log << Messages::endl;
}


void
Messages::close_log_file(void)
{
  time_t now;
  time(&now);
  _log << Messages::endl;
  _log << "Log end: " << ctime(&now) << Messages::endl;
  _log << Messages::endl;
  _log.close();
}



void
Messages::warning(const string& msg)
{
  TeeStream ts(cout, _log);
  ts << Messages::endl;
#ifdef _WIN32
  HANDLE hstdout = GetStdHandle(STD_OUTPUT_HANDLE);
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  GetConsoleScreenBufferInfo(hstdout, &csbi);
  SetConsoleTextAttribute(hstdout, FOREGROUND_GREEN | FOREGROUND_BLUE);
#else
  cout << yellowb;
#endif
  _log << "*** ";
  ts << "Warning: ";
#ifdef _WIN32
  SetConsoleTextAttribute(hstdout, csbi.wAttributes);
#else
  cout << normal;
#endif
  ts << msg << endl << flush;

  if (interactive())
  {
    cout << "press any key to continue ...";
    cin.get();
  }

  _warning_count++;
}



void
Messages::error(const string& msg)
{
  TeeStream ts(cerr, _log);
  ts << Messages::endl;
#ifdef _WIN32
  HANDLE hstdout = GetStdHandle(STD_OUTPUT_HANDLE);
  CONSOLE_SCREEN_BUFFER_INFO csbi;
  GetConsoleScreenBufferInfo(hstdout, &csbi);
  SetConsoleTextAttribute(hstdout, FOREGROUND_RED);
#else
  cerr << redb;
#endif
  _log << "*** ";
  ts << "ERROR: ";
#ifdef _WIN32
  SetConsoleTextAttribute(hstdout, csbi.wAttributes);
#else
  cout << normal;
#endif
  ts << msg << endl << flush;

  if (interactive())
  {
    cout << "press any key to continue ...";
    cin.get();
  }

  _error_count++;
}



void
Messages::info(const string& msg, bool newline)
{
  static bool contd = false;

  vector<string> lines;
  Utils::tokenize(msg, lines, "\n");

  TeeStream ts(cout, _log);

  size_t nl = lines.size();

  for (size_t l = 0; l < nl; l++)
  {
    if (!contd)
      for (int i = 0; i < _indent * _indent_width; i++)
        ts << " ";
    ts << lines[l];
    if (newline || (l < nl - 1)) ts << endl;

    ts << flush;
  }

  if (newline) contd = false;
  else contd = true;

}


void
Messages::debug(const string& msg)
{
#ifdef DEBUG
  TeeStream ts(cerr, _log);
  ts << "DEBUG: " << msg <<  endl << flush;
#else
  ignore_unused_variable(msg);
#endif
}


void
Messages::newline(void)
{
  TeeStream ts(cout, _log);
  ts << endl << flush;
}

void
Messages::frameline(const std::string& start, const char ch, const std::string& name)
{
  //ostringstream os;
  //os << start;
  //os.fill(ch);
  //os.width(Messages::available_width() - start.length());
  //os << "" << Messages::endl;

  int w = Messages::available_width();
  int wn = name.length();
  int wc = wn > 0 ? 2 : 0;
  string s;
  s.reserve(w);
  s += start;
  w -= start.length() + 2*wc + name.length();
  for (int i = 0; i < w / 2; i++)
    s += ch;
  if (wn > 0) s += " (";
  s += name;
  if (wn > 0) s += ") ";
  for (int i = 0; i < (w / 2 + w % 2); i++)
    s += ch;
  Messages::info(s);

}

