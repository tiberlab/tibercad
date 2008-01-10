// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "Control.h"
#include "DLLoader.h"
#include "TiberCad.h"


#include <iostream>



using namespace std;


// Will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  cout << "TiberCAD version 0.1.0" << endl << endl;

  // take input file from command line or ask for it
  string inputfile;
  if (argc > 1)
    inputfile = string(argv[1]);
  else
  {
    //cerr << "Usage: tibercad <inputfile>" << endl;
    //return 1;

    cout << "input file: ";
    cin >> inputfile;
    boost::algorithm::trim(inputfile);
    cout << endl;
  }

  // Set up some path
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
    cout << "ERROR: TiberCAd crashed for unknown reason." << endl;
  }

  return TiberCad::cleanup();
}

