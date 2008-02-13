// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "tiber_config.h"

#include "Control.h"
#include "DLLoader.h"
#include "TiberCad.h"
#include "License.h"


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
    cerr << "Usage: tibercad <inputfile>" << endl << endl;
# ifdef CYGWIN
    cout << "press any key ...";
    cin.get();
    return 1;
# endif
#endif

  }

  // do some preparation
  {
    // the TiberCAD root
    string tiberroot;
    char* root = getenv("TIBERCADROOT");
    if (root != NULL)
      tiberroot = string(root);

    // Set up search path for libraries
    DLLoader::set_library_path(tiberroot + "/lib/tibermodels");

    char* modelpath = getenv("TIBERMODELPATH");
    if (modelpath != NULL)
      DLLoader::prepend_to_library_path(modelpath);

    //DLLoader::prepend_to_library_path(".");

    // Set up search path for materials

#ifdef LICENSE_CHECK
    // check the license
    bool lic_ok = false;
    string default_licfile("tibercad.lic");
    char* licfile = getenv("TIBERLICENSEFILE");
    if (licfile != NULL)
      lic_ok = License::check_license(licfile);

    if (!lic_ok)
      lic_ok = License::check_license(default_licfile);

    if (!lic_ok)
      lic_ok = License::check_license(tiberroot + "/" + default_licfile);

    if (!lic_ok)
      lic_ok = License::check_license(tiberroot + "/license/" + default_licfile);

    if (!lic_ok)
    {
      cerr << "Sorry, cannot start TiberCAD as I could not find a valid "
        << "license." << endl;
# ifdef CYGWIN
      cout << endl << "press any key ...";
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
#ifdef CYGWIN
    cout << "press any key ...";
    cin.get();
#endif
  }
  catch (...)
  {
    cout << "ERROR: TiberCAD crashed for unknown reason." << endl;
#ifdef CYGWIN
    cout << "press any key ...";
    cin.get();
#endif
  }

  return TiberCad::cleanup();
}

