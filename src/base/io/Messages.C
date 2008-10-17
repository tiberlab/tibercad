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

  string _warning = string(yellow) + "Warning: " + string(white);
  string _error   = string(red) + "ERROR: " + string(white);
  string _reset   = normal;

}

void
Messages::warning(const string& msg)
{
  cout << _warning << msg << _reset << endl;
}



void
Messages::error(const string& msg)
{
  cerr << _error << msg << _reset << endl;
}



void
Messages::info(const string& msg)
{
  cerr << msg << endl;
}


