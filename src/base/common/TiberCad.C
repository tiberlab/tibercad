// $Id$

#include "tiber_config.h"
#include "tiber_version.h"
#include "svnrevision.h"
#include "TiberCad.h"
#include "Control.h"
#include "EigenSolver.h"
#include "Database.h"
#include "DLLoader.h"
#include "Utils.h"
#include "Messages.h"
#include "InitFailedException.h"

#include "libmesh.h"
#include "petscerror.h"


#ifndef ARCH
#error "Architecture has to be specified on the command line as string"
#endif

// we hand an empty command line to underlying libraries
namespace
{
  int __empty_argc = 1;
  char** __empty_argv = new char*[1];
  char __executable[] = "tibercad";
}


std::list<std::string>
TiberCad::_filename_suffix;


Control*
TiberCad::_control = NULL;


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




TiberCad::TiberCad(void) :
  _libmeshinit(NULL)
{
  if (_object_counter > 0)
    throw InitFailedException("Only one TiberCAD instance may exist in a process!");

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



std::string
TiberCad::arch_string(void)
{
  return std::string(ARCH);
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
TiberCad::init(const std::string& inputfile)
{
  // read TIBERCADROOT from environment
  char* root = getenv("TIBERCADROOT");
  if (root != NULL)
    _tiberroot = std::string(root);

  Messages::debug("Using TIBERCADROOT=" + std::string(root));

  if (_tiberroot.size() != 0)
  {
    // setup default database search path
    Database::set_default_search_path(_tiberroot + "/materials");

    // setup DLLoader paths
#ifdef _WIN32
    DLLoader::set_library_path(_tiberroot + "/modules");
#elif defined(__APPLE__)
    DLLoader::set_library_path(_tiberroot + "/modules");
    DLLoader::set_library_path(_tiberroot + "/" + ARCH + "/modules");
#else
    DLLoader::append_to_library_path(_tiberroot + "/" + ARCH + "/modules");
#ifdef DEBUG
    DLLoader::prepend_to_library_path(_tiberroot + "/" + ARCH + "/modules/debug");
#endif
#endif
    char* modelpath = getenv("TIBERMODULEPATH");
    if (modelpath != NULL)
      DLLoader::prepend_to_library_path(modelpath);
  }


  // to the libraries we hand empty cmdline!
  __empty_argv[0] = __executable;


  MPI_Comm local_comm;
  int proc_id, p;

  /* starts MPI */
  MPI_Init(&__empty_argc, &__empty_argv);
  /* get current process id */
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);    /* get number of processes */
  MPI_Comm_split(MPI_COMM_WORLD, proc_id, 0, &local_comm);

  int local_id;
  MPI_Comm_rank(local_comm, &local_id);
  MPI_Barrier(MPI_COMM_WORLD);

  //if (proc_id == 1) sleep(2);
  //std::cerr << "global: rank = " << proc_id << " np = " << p << std::endl;
  //std::cerr << "local : rank = " << local_id << std::endl;

  // prepare libMesh
  _libmeshinit = new LibMeshInit(__empty_argc, __empty_argv, local_comm);

  // prepare EigenSolver
  EigenSolver::slepc_init(__empty_argc, __empty_argv);//, local_comm);

  PetscPopSignalHandler();


  // now create a TiberCAD Control object
  _control = new Control();
 
  std::string infile(inputfile);
#if defined(__CYGWIN__)
    // we first convert the filename to something more UNIX like
    Utils::convert_win32_path_to_posix(infile);
#endif

  _control->set_inputfile(infile);
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

  MPI_Finalize();

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
