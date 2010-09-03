// $Id$

#include "tiber_config.h"
#include "tiber_version.h"
#include "svnrevision.h"
#include "TiberCad.h"
#include "Control.h"
#include "EigenSolver.h"
#include "Database.h"
#include "DLLoader.h"

#include "libmesh.h"

namespace
{
  Control* _control = NULL;

  LibMeshInit* _libmeshinit;
}

char**
TiberCad::_cmdline_argv = 0;

int
TiberCad::_cmdline_argc = 0;


std::string
TiberCad::_tiberroot = "";


const int
TiberCad::_MajorVersion = TIBERMAJORVERSION;


const int
TiberCad::_MinorVersion = TIBERMINORVERSION;


const int
TiberCad::_SubMinorVersion = TIBERSUBMINORVERSION;


const int
TiberCad::_SvnRevision = SVNREVISION;



std::string
TiberCad::version_string(bool include_svn_release)
{
  std::ostringstream os;
  os << _MajorVersion << "." << _MinorVersion << "." << _SubMinorVersion;
  if (include_svn_release)
    os << " rev. " << _SvnRevision;

  return os.str();
}

int
TiberCad::major_version(void)
{
  return _MajorVersion;
}


int
TiberCad::minor_version(void)
{
  return _MinorVersion;
}



int
TiberCad::subminor_version(void)
{
  return _SubMinorVersion;
}



int
TiberCad::software_revision(void)
{
  return _SvnRevision;
}




void
TiberCad::init(int argc, char** argv)
{
  // perhaps we can use the commandline arguments somewhere?
  _cmdline_argc = argc;
  _cmdline_argv = argv;

  // read TIBERCADROOT from environment
  char* root = getenv("TIBERCADROOT");
  if (root != NULL)
    _tiberroot = std::string(root);

  if (_tiberroot.size() != 0)
  {
    // setup default database search path
    Database::set_default_search_path(_tiberroot + "/materials");

    // setup DLLoader paths
    DLLoader::set_library_path(_tiberroot + "/lib/modules");
#ifdef DEBUG
    DLLoader::prepend_to_library_path(_tiberroot + "/lib/debug/modules");
#endif
    char* modelpath = getenv("TIBERMODELPATH");
    if (modelpath != NULL)
      DLLoader::prepend_to_library_path(modelpath);
  }


  // prepare libMesh
  //libMesh::init(cmdline_argc, cmdline_argv);
  _libmeshinit = new LibMeshInit(_cmdline_argc, _cmdline_argv);

  // prepare EigenSolver
  EigenSolver::slepc_init(_cmdline_argc, _cmdline_argv);


  // now create a TiberCAD Control object
  _control = new Control();
}


int
TiberCad::cleanup(void)
{
  delete _control;

  // close EigenSolver
  EigenSolver::slepc_done();

  // close libMesh and return
  delete _libmeshinit;
  //return libMesh::close();
  return 0;
}


Control&
TiberCad::get_control(void)
{
  assert(_control != NULL);
  return *_control;
}


const std::string&
TiberCad::get_output_format(void)
{
  return get_control().get_output_format();
}


const std::string&
TiberCad::get_output_dir(void)
{
  return get_control().get_output_dir();
}



const std::string&
TiberCad::get_filename_suffix(void)
{
  return get_control().get_filename_suffix();
}


void
TiberCad::clear_filename_suffix(void)
{
  get_control().clear_filename_suffix();
}


void
TiberCad::append_to_filename_suffix(const std::string& suffix)
{
  get_control().append_to_filename_suffix(suffix);
}



void
TiberCad::prepend_to_filename_suffix(const std::string& suffix)
{
  get_control().prepend_to_filename_suffix(suffix);
}


void
TiberCad::drop_first_filename_suffix(void)
{
  get_control().drop_first_filename_suffix();
}


void
TiberCad::drop_last_filename_suffix(void)
{
  get_control().drop_last_filename_suffix();
}
