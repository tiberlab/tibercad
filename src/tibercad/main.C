// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "tiber_config.h"

#include "Control.h"
#include "DLLoader.h"
#include "TiberCad.h"
#include "License.h"
#include "Utils.h"


#include <iostream>
#ifdef HAVE_LIBREADLINE
# ifdef HAVE_READLINE_READLINE_H
#  include <readline/readline.h>
# else
#  include <readline.h>
# endif
#endif



using namespace std;


// Will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  cout << "TiberCAD version " << TIBERVERSION << endl << endl;


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
    cout << "Usage:" << endl
      << "  from command line: tibercad <inputfile>" << endl
      << "  or double click on inputfile" << endl << endl;
    cout << "press Enter ...";
    cin.get();
# else
    cout << "Usage: tibercad <inputfile>" << endl << endl;
# endif
    return 1;
#endif

  }

  // do some preparation
  {
#ifdef CYGWIN
    // in windows argv[0] is the absolute path
    char* root = getenv("TIBERCADROOT");
    if (root == NULL)
    {
      string exepath(Utils::dirname(argv[0]));
      setenv("TIBERCADROOT", exepath.c_str(), 1);
    }
#endif
    // the TiberCAD root
    //char* root = getenv("TIBERCADROOT");
    //if (root != NULL)
    //  TiberCAD::tiberroot = string(root);


    // Set up search path for libraries
    //DLLoader::set_library_path(tiberroot + "/lib/tibermodels");

    //char* modelpath = getenv("TIBERMODELPATH");
    //if (modelpath != NULL)
    //  DLLoader::prepend_to_library_path(modelpath);

    //DLLoader::prepend_to_library_path(".");

    // Set up search path for materials

#ifdef LICENSE_CHECK
    // check the license
    if (!License::check_license())
    {
      cerr << "Sorry, cannot start TiberCAD as I could not find a valid "
        << "license." << endl;
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
    cout << "Simulation finished..." << endl << "Goodbye" << endl;

  }
  catch (exception& e)
  {
    cout << "ERROR: " << e.what() << endl;
  }
  catch (...)
  {
    cout << "ERROR: TiberCAD crashed for unknown reason." << endl;
  }
#ifdef CYGWIN
  cout << "press Enter ...";
  cin.get();
#endif

  return TiberCad::cleanup();
}

