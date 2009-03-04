// $Id$

#include "Messages.h"
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


string
Messages::_warning = string(yellow) + "Warning: " + string(normal);

string
Messages::_error   = string(red) + "ERROR: " + string(normal);

string
Messages::_reset   = normal;

string
Messages::_endl    = string(normal) + "\n";




void
Messages::warning(const string& msg)
{
  cout << _warning << msg << _reset << _endl << flush;
}



void
Messages::error(const string& msg)
{
  cerr << _error << msg << _reset << _endl << flush;
}



void
Messages::info(const string& msg)
{
  vector<string> lines;
  Utils::tokenize(msg, lines, "\n");

  for (int l = 0; l < lines.size(); l++)
  {
    for (int i = 0; i < _indent; i++)
      cout << "  ";
    cout << lines[l] << _endl << flush;
  }
}


void
Messages::debug(const std::string& msg)
{
#ifdef DEBUG
  cerr << "DEBUG: " << msg <<  endl << flush;
#else
  ignore_unused_variable(msg);
#endif
}


void
Messages::newline(void)
{
  cout << endl << flush;
}
