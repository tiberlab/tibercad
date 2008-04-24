// $Id$

#include "tiber_config.h"
#include "TiberCad.h"
#include "EigenSolver.h"

#include "libmesh.h"


char**
TiberCad::cmdline_argv = 0;

int
TiberCad::cmdline_argc = 0;


std::string
TiberCad::tiberroot = "";


const int
TiberCad::MajorVersion = 0;


const int
TiberCad::MinorVersion = 1;


const int
TiberCad::SubMinorVersion = 0;


const int
TiberCad::SvnRevision = SVNREVISION;


std::string
TiberCad::TiberCadVersion(bool include_svn_release)
{
  std::ostringstream os;
  os << MajorVersion << "." << MinorVersion << "." << SubMinorVersion;
  if (include_svn_release)
    os << "." << SvnRevision;

  return os.str();
}



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
