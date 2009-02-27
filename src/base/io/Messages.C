// $Id$

#include "Messages.h"


using namespace std;


namespace
{
  const char* red    = "\033[1;31m";
  const char* yellow = "\033[1;33m";
  const char* blue   = "\033[1;34m";
  const char* white  = "\033[1;37m";
  const char* normal = "\033[0m";
}


string
Messages::_warning = string(yellow) + "Warning: " + string(normal);

string
Messages::_error   = string(red) + "ERROR  : " + string(normal);

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
  cerr << msg << _endl << flush;
}


