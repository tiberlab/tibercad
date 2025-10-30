// $Id$

#include "tiber_config.h"
#include "TiberCad.h"
#include "control/Control.h"
#include "EigenSolver.h"
#include "Database.h"
#include "DLLoader.h"
#include "Utils.h"
#include "Messages.h"
#include "InitFailedException.h"

#include "libmesh/libmesh.h"
#include "petscsys.h"
#include "petscerror.h"

//#ifdef _WIN32
 #include "omp.h"
//#endif


#ifndef ARCH
#error "Architecture has to be specified on the command line as string"
#endif

// we hand an empty command line to underlying libraries
namespace
{
  int __empty_argc = 2;
  char** __empty_argv = new char*[2];
  char __executable[] = "tibercad";
  char __option[] = "--node-major-dofs";
}


std::list<std::string>
TiberCad::_filename_suffix;


Control*
TiberCad::_control = NULL;


unsigned int
TiberCad::_object_counter = 0;

bool
TiberCad::_is_initialized = false;

std::string
TiberCad::_tiberroot = "";


const std::string
TiberCad::_git_revision = TC_REVISION;

const std::string
TiberCad::_git_modified = TC_MODIFIED;

const std::string
TiberCad::_compilation_date = TC_COMPDATE;

const std::string
TiberCad::_compilation_system = TC_COMPSYS;


libMesh::Parallel::Communicator
TiberCad::_mpi_comm;


TiberCad::TiberCad(MPI_Comm mpi_comm) :
  _libmeshinit(NULL)
{
  if (_object_counter > 0)
    throw InitFailedException("Only one TiberCAD instance may exist in a process!");

  _object_counter++;

  _mpi_comm.duplicate(mpi_comm);
}


TiberCad::~TiberCad(void)
{
  cleanup();
}



libMesh::Parallel::Communicator&
TiberCad::get_mpi_comm()
{
  return(_mpi_comm);
}


std::string
TiberCad::version_string(bool include_compilation_date)
{
  std::ostringstream os;
  os << _git_revision;
  if (include_compilation_date)
    os << " compiled on " << _compilation_date;

  return os.str();
}

const std::string&
TiberCad::compilation_date(void)
{
  return _compilation_date;
}


const std::string&
TiberCad::compilation_system(void)
{
  return _compilation_system;
}


const std::string&
TiberCad::last_modification(void)
{
  return _git_modified;
}



std::string
TiberCad::arch_string(void)
{
  return std::string(ARCH);
}



const std::string&
TiberCad::software_revision(void)
{
  return _git_revision;
}



void
TiberCad::init(const std::string& inputfile)
{
  if (_is_initialized)
    return;

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
  __empty_argv[1] = __option;

  // Note: it is the responsability of the user to define OMP_NUM_THREADS correctly
  // in case of MPI parallel runs
  int omp_procs = omp_get_max_threads();

  // prepare libMesh
  _libmeshinit = new libMesh::LibMeshInit(__empty_argc, __empty_argv, _mpi_comm.get(), omp_procs);

  // do not pass all arguments to slepc
  __empty_argc = 1;
  EigenSolver::slepc_init(__empty_argc, __empty_argv, _mpi_comm.get());

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

  _is_initialized = true;
}



void
TiberCad::run(void)
{
  if (!_is_initialized)
    throw InitFailedException("TiberCAD library is not initialized!");

  _control->run_simulation();
}



void
TiberCad::cleanup(void)
{
  if (!_is_initialized)
    return;

  // close EigenSolver
  EigenSolver::slepc_done();

  //delete _control;

  // close libMesh and return
  delete _libmeshinit;

  _mpi_comm.clear();

  delete [] __empty_argv;

  _is_initialized = false;
}


double
TiberCad::get_global_time(void)
{
  return _control->get_time();
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
