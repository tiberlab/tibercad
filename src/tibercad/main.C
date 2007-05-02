// $Id$


#include "Control.h"
#include "DLLoader.h"

#include "libmesh.h"

#include <iostream>

#include "EigenSolver.h"


using namespace std;


// This is the first TiberCAD main !!!
// It will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  string inputfile;
  if (argc > 1)
    inputfile = string(argv[1]);
  else
  {
    cerr << "Usage: tibercad <inputfile>" << endl;
    return 1;
  }

  cout << "TiberCAD version 0.1.0" << endl << endl;

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

  libMesh::init(argc, argv);
  {
    EigenSolver::slepc_init();

    Control control(inputfile);
    
    control.init();
    control.run_simulation();

    EigenSolver::slepc_done();
  }

  cout << "Simulation finished..." << endl << "Goodbye" << endl;

  return libMesh::close();
}

