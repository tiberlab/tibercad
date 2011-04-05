// $Id$

#include "tiber_config.h"
#include "tiber_version.h"
#include "svnrevision.h"
#include "TiberCad.h"
#include "Control.h"
#include "EigenSolver.h"
#include "Database.h"
#include "DLLoader.h"
#include "InitFailedException.h"

#include "libmesh.h"


#ifndef stringify
#define stringify(a) #a
#endif
#ifndef xstr
#define xstr(a) stringify(a)
#endif


// we hand an empty command line to underlying libraries
namespace
{
  int __empty_argc = 1;
  char** __empty_argv = new char*[1];
}


std::list<std::string>
TiberCad::_filename_suffix;


Control*
TiberCad::_control = NULL;


char**
TiberCad::_cmdline_argv = 0;

int
TiberCad::_cmdline_argc = 0;

unsigned int
TiberCad::_object_counter = 0;


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




TiberCad::TiberCad(int argc, char** argv) :
  _libmeshinit(NULL)
{
  if (_object_counter > 0)
    throw InitFailedException("Only one TiberCAD instance may exist in a process!");

  _cmdline_argc = argc;
  _cmdline_argv = argv;

  _object_counter++;
}


TiberCad::~TiberCad(void)
{
  cleanup();
}



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
TiberCad::init(void)
{
  // read TIBERCADROOT from environment
  char* root = getenv("TIBERCADROOT");
  if (root != NULL)
    _tiberroot = std::string(root);

  if (_tiberroot.size() != 0)
  {
    // setup default database search path
    Database::set_default_search_path(_tiberroot + "/materials");

    // setup DLLoader paths
#ifdef CYGWIN
    DLLoader::set_library_path(_tiberroot + "/lib/modules");
#else
    DLLoader::set_library_path(_tiberroot + "/" + xstr(ARCH) + "/lib/modules");
#ifdef DEBUG
    DLLoader::prepend_to_library_path(_tiberroot + "/" + xstr(ARCH) + "/lib/debug/modules");
#endif
#endif
    char* modelpath = getenv("TIBERMODULEPATH");
    if (modelpath != NULL)
      DLLoader::prepend_to_library_path(modelpath);
  }


  // to the libraries we hand empty cmdline!
  __empty_argv[0] = _cmdline_argv[0];


  // prepare libMesh
  _libmeshinit = new LibMeshInit(__empty_argc, __empty_argv);

  // prepare EigenSolver
  EigenSolver::slepc_init(__empty_argc, __empty_argv);


  // now create a TiberCAD Control object
  _control = new Control();
 
  std::string inputfile(_cmdline_argv[1]);
#ifdef CYGWIN
    // we first convert the filename to something more UNIX like
    Utils::convert_win32_path_to_posix(inputfile);
#endif

  _control->set_inputfile(inputfile);
  _control->init();
}



void
TiberCad::run(void)
{
  _control->run_simulation();
}



void
TiberCad::cleanup(void)
{
  delete _control;

  // close EigenSolver
  EigenSolver::slepc_done();

  // close libMesh and return
  delete _libmeshinit;

  delete [] __empty_argv;
}




const std::string&
TiberCad::get_output_format(void)
{
  return _control->get_output_format();
}


const std::string&
TiberCad::get_output_dir(void)
{
  return _control->get_output_dir();
}



std::string
TiberCad::get_filename_suffix(void)
{
  std::string suffix_str;
  std::list<std::string>::const_iterator it(_filename_suffix.begin());
  const std::list<std::string>::const_iterator end(_filename_suffix.end());
  for ( ; it != end; ++it)
    suffix_str += "_" + *it;

  return suffix_str;
}


void
TiberCad::clear_filename_suffix(void)
{
  _filename_suffix.clear();
}


void
TiberCad::append_to_filename_suffix(const std::string& suffix)
{
  if (suffix.size() != 0)
    _filename_suffix.push_back(suffix);
}



void
TiberCad::prepend_to_filename_suffix(const std::string& suffix)
{
  if (suffix.size() != 0)
    _filename_suffix.push_front(suffix);
}


void
TiberCad::drop_first_filename_suffix(void)
{
  if (_filename_suffix.size() != 0)
    _filename_suffix.pop_front();
}


void
TiberCad::drop_last_filename_suffix(void)
{
  if (_filename_suffix.size() != 0)
    _filename_suffix.pop_front();
}
