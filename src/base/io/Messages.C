// $Id$

#include "Messages.h"
#include "TeeStream.h"
#include "TypeDefs.h"
#include "Utils.h"

#include <vector>


namespace
{
  const char* red     = "\033[1;31m";
  const char* yellow  = "\033[1;33m";
  const char* redb    = "\033[1;41m";
  const char* yellowb = "\033[1;43m";
  const char* blue    = "\033[1;34m";
  const char* white   = "\033[1;37m";
  const char* normal  = "\033[0m";
}


const std::string
Messages::endl = "\n";

const int
Messages::_max_width = 80;

int
Messages::_indent = 0;


std::ofstream
Messages::_log;


void
Messages::set_log_file(const std::string& logfile)
{
  _log.open(logfile.c_str());
}


void
Messages::close_log_file(void)
{
  _log.close();
}



void
Messages::warning(const std::string& msg)
{
  TeeStream ts(std::cout, _log);
  ts << Messages::endl;
  std::cout << yellowb;
  _log << "*** ";
  ts << "Warning:";
  std::cout << normal << " ";
  ts << msg << endl << std::flush;
}



void
Messages::error(const std::string& msg)
{
  TeeStream ts(std::cerr, _log);
  ts << Messages::endl;
  std::cerr << redb;
  _log << "*** ";
  ts << "ERROR:";
  std::cerr << normal << " ";
  ts << msg << endl << std::flush;
}



void
Messages::info(const std::string& msg)
{
  std::vector<std::string> lines;
  Utils::tokenize(msg, lines, "\n");

  TeeStream ts(std::cout, _log);

  for (size_t l = 0; l < lines.size(); l++)
  {
    for (int i = 0; i < _indent; i++)
      ts << "  ";
    ts << lines[l] << endl << std::flush;
  }
}


void
Messages::debug(const std::string& msg)
{
#ifdef DEBUG
  TeeStream ts(std::cerr, _log);
  ts << "DEBUG: " << msg <<  endl << std::flush;
#else
  ignore_unused_variable(msg);
#endif
}


void
Messages::newline(void)
{
  TeeStream ts(std::cout, _log);
  ts << endl << std::flush;
}
