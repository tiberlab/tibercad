// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "Control.h"
#include "DLLoader.h"
#include "TiberCad.h"
#include "License.h"
#include "Utils.h"
#include "Messages.h"

#include "tiber_config.h"


#include <iostream>
#ifdef HAVE_LIBREADLINE
# ifdef HAVE_READLINE_READLINE_H
#  include <readline/readline.h>
# else
#  include <readline.h>
# endif
#endif

#ifdef CYGWIN
# include <windows.h>
#endif


using namespace std;


// Will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  // take input file from command line or ask for it
  string inputfile;
  if (argc > 1)
    inputfile = string(argv[1]);
  else
  {

#ifdef HAVE_LIBREADLINE
    char *line = readline ("Enter input file: ");
    inputfile = string(line);
    free(line);
    boost::algorithm::trim(inputfile);
    cout << endl;
#else
# ifdef CYGWIN
    cout << endl << "Usage:" << endl
      << "  from command line: tibercad <inputfile>" << endl
      << "  or double click on inputfile" << endl << endl;
    cout << "press Enter ...";
    cin.get();
# else
    cout << endl << "Usage: tibercad <inputfile>" << endl << endl;
# endif
    return 1;
#endif

  }

  Messages::set_log_file(Utils::basename(inputfile) + ".log");

  Messages::newline();
  Messages::info("TiberCAD version " + TiberCad::TiberCadVersion());
  Messages::newline();

  Messages::info("Input file: " + inputfile);
  Messages::newline();

  // do some preparation
  {
#ifdef CYGWIN
    // we first convert the filename to something more UNIX like
    Utils::convert_path_to_unix(inputfile);

    // in windows argv[0] is the absolute path
    char* root = getenv("TIBERCADROOT");
    if (root == NULL)
    {
      const size_t bufsize = 1024;
      char buffer[bufsize];
      if (!GetModuleFileName(NULL, buffer, bufsize))
        Messages::warning("Problems detecting installation path.");
      string program(buffer);
      Utils::convert_path_to_unix(program);
      string exepath(Utils::dirname(program));
      setenv("TIBERCADROOT", exepath.c_str(), 1);
    }
#endif


#ifdef LICENSE_CHECK
    // check the license
    if (!License::check_license())
    {
      Messages::info("Sorry, cannot start TiberCAD as I could not find "
          "a valid license.");
# ifdef CYGWIN
      cout << endl << "press Enter ...";
      cin.get();
# endif
      return 1;
    }
#endif

  }

  //
  // here begins real TiberCAD
  //
  TiberCad::init(argc, argv);
  try {

    Control control(inputfile);
    
    control.init();
    control.run_simulation();

    Messages::info("Simulation finished...");
    Messages::info("Goodbye");

  }
  catch (exception& e)
  {
    Messages::error(e.what());
  }
  catch (...)
  {
    Messages::error("TiberCAD crashed for unknown reason.");
  }
#ifdef CYGWIN
  cout << "press Enter ...";
  cin.get();
#endif

  return TiberCad::cleanup();
}

