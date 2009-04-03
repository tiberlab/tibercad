// $Id$

#include "Messages.h"
#include "TeeStream.h"
#include "TypeDefs.h"
#include "Utils.h"

#include <vector>

using namespace std;


namespace
{
  const char* red    = "\033[1;31m";
  const char* yellow = "\033[1;33m";
  const char* blue   = "\033[1;34m";
  const char* white  = "\033[1;37m";
  const char* normal = "\033[0m";
}


int
Messages::_indent = 0;


ofstream
Messages::_log;


void
Messages::set_log_file(const string& logfile)
{
  _log.open(logfile.c_str());
}


void
Messages::close_log_file(void)
{
  _log.close();
}



void
Messages::warning(const string& msg)
{
  TeeStream ts(cout, _log);
  cout << yellow;
  _log << "*** ";
  ts << "Warning: ";
  cout << normal;
  ts << msg << endl << flush;
}



void
Messages::error(const string& msg)
{
  TeeStream ts(cerr, _log);
  cerr << red;
  _log << "*** ";
  ts << "ERROR: ";
  cerr << normal;
  ts << msg << endl << flush;
}



void
Messages::info(const string& msg)
{
  vector<string> lines;
  Utils::tokenize(msg, lines, "\n");

  TeeStream ts(cout, _log);

  for (int l = 0; l < lines.size(); l++)
  {
    for (int i = 0; i < _indent; i++)
      ts << "  ";
    ts << lines[l] << endl << flush;
  }
}


void
Messages::debug(const std::string& msg)
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
