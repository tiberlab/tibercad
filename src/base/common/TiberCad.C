// $Id$

#include "TiberCad.h"
#include "EigenSolver.h"

#include "libmesh.h"


char**
TiberCad::cmdline_argv = 0;

int
TiberCad::cmdline_argc = 0;


std::string
TiberCad::tiberroot = "";


void
TiberCad::init(int argc, char** argv)
{
  // perhaps we can use the commandline arguments somewhere?
  cmdline_argc = argc;
  cmdline_argv = argv;

  // read TIBERCADROOT from environment
  char* root = getenv("TIBERCADROOT");
  if (root != NULL)
    tiberroot = std::string(root);

  
  // prepare libMesh
  libMesh::init(cmdline_argc, cmdline_argv);
  
  // prepare EigenSolver
  EigenSolver::slepc_init(cmdline_argc, cmdline_argv);
}


int
TiberCad::cleanup(void)
{

  // close EigenSolver
  EigenSolver::slepc_done();

  // close libMesh and return
  return libMesh::close();
}
