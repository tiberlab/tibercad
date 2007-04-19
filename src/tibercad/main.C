// $Id$


#include "Control.h"

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

  libMesh::init(argc, argv);
  {
    slepc_init();

    Control control(inputfile);
    
    control.init();
    control.run_simulation();

    slepc_done();
  }

  cout << "Simulation finished..." << endl << "Goodbye" << endl;

  return libMesh::close();
}

