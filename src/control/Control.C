// $Id$


#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "base/io/InputParser.h"
#include "Control.h"
#include "Database.h"
#include "base/io/DLLoader.h"
#include "Utils.h"
#include "Variable.h"
#include "SignalGenerator.h"
#include "Messages.h"
#include "SimulationOptions.h"
#include "Device.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "PhysicalModel.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationInterface.h"
#include "RuntimeException.h"
#include "SolveFailedException.h"
#include "ModelErrorException.h"
#include "InitFailedException.h"

#include <sstream>
#include <vector>
#include <set>

#include <csignal>


using namespace std;



/*
Control*
Control::SignalHandler::_ctrl = NULL;

Control::SignalHandler::SigAction
Control::SignalHandler::_old_int_action;


void
Control::SignalHandler::set_control(Control* ctrl)
{
  _ctrl = ctrl;
}


void
Control::SignalHandler::activate_sigint(void)
{
  sigaction(SIGINT, NULL, &_old_int_action);

  SigAction new_action;
  new_action.sa_flags = SA_RESTART;
  sigemptyset(&new_action.sa_mask);
  new_action.sa_handler = sigint;

  if (_old_int_action.sa_handler != SIG_IGN)
    sigaction(SIGINT, &new_action, NULL);
}




void
Control::SignalHandler::deactivate_sigint(void)
{
  sigaction(SIGINT, &_old_int_action, NULL);
}



void
Control::SignalHandler::sigint(int sig)
{
  // A second Ctrl-C should silently quit
  deactivate_sigint();
  sigset_t block_int;
  sigemptyset(&block_int);
  sigaddset(&block_int, SIGINT);
  sigprocmask(SIG_UNBLOCK, &block_int, NULL);

  cerr << "\nDo you really want to quit TiberCAD [n]? ";
  char buffer[16];
  cin.get(buffer, 16, '\n');
  string ans(buffer);
  cin.clear();
  cin.ignore(256, '\n');

  boost::algorithm::trim(ans);
  boost::algorithm::to_lower(ans);

  if ((ans == "y") || (ans == "yes"))
  {
    raise(sig);
  }

  // re-block the signal
  sigprocmask(SIG_BLOCK, &block_int, NULL);

  // reactivate the handler
  activate_sigint();
}
*/


Control::Control(void)
{

}



Control::~Control(void)
{
  SimulationInterface::SimulationIterator
    simit(SimulationInterface::simulations_begin());
  const SimulationInterface::SimulationIterator
    simend(SimulationInterface::simulations_end());

  while (simit != simend)
  {
    SimulationInterface* sim = *simit;
    ++simit;
    SimulationInterface::destroy(sim);
  }

  // clear all variables
  VariableValue::clear_all();

  // This in the past gave a problem with some MPI comm
  Device::destroy(_device);

}





void
Control::init(void)
{

  // TODO this is currently broken
  // setup the signal handler
  //SignalHandler::set_control(this);

  // we want to intercept SIGINT (Ctrl-C)
  //SignalHandler::activate_sigint();


  {
    // get the global rank
    // we then define MPI_PROC to the global rank so we can
    // use this in the input file to change e.g. output directory
    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    ostringstream os;
    os << proc_id;
    InputParser::add_defined("MPI_PROC", os.str());
  }

  Messages::set_communicator(TiberCad::get_mpi_comm(), 0);

  // I would like to define preprocessor variables from device
  // MPI comm groups, but for that the input parser is not flexible enough
  // Therefore at the moment we can do a trick and parse twice
  // First parsing, and create a fake device:
  ModelOptions input;
  InputParser parser;
  parser.parse_file(_inputfile, input);

  ModelOptions::submodel_iterator it = input.submodels_begin("Device");
  if (it == input.submodels_end("Device"))
    throw InitFailedException("\'Device\' block missing in input file.");


  // this will be the communicator passed to Device constructor
  libMesh::Parallel::Communicator dev_comm;

  // here we set up the communicator for the device
  {
    ModelOptions& dev_opts = it->second;
    libMesh::Parallel::Communicator &comm = TiberCad::get_mpi_comm();
    unsigned int nodes_per_device = comm.size();

    int color = 0;

    if (dev_opts.has_submodel("Parallel"))
    {
      const ModelOptions &mpi_opts = dev_opts.submodels_begin("Parallel")->second;

      // nodes_per_device is the number of processes to be used for calculation of a device.
      // If e.g. comm.size() = 4 and nodes_per_device = 2, then simulation on two independent
      // devices is performed, each parallelized on two processes
      nodes_per_device = mpi_opts.get_option("mpi_processes_per_device", nodes_per_device);
    }

    if (nodes_per_device < comm.size())
    {
      int proc_id = comm.rank();
      color = proc_id / nodes_per_device;

      comm.split(color, 0, dev_comm);
    }
    else
    {
      // TODO not sure if we should duplicate here?
      //_mpi_comm = TiberCad::get_mpi_comm();
      dev_comm.duplicate(TiberCad::get_mpi_comm());

      if (nodes_per_device > comm.size())
        throw InitFailedException("Too many MPI nodes requested for device");
    }

    // set a special variable to get fine control from input file
    // e.g. one might use MPI_DEV_KEY to change output directories per device
    ostringstream os;
    os << color;
    InputParser::add_defined("MPI_DEV_KEY", os.str(), false);
  }

  // second parsing, now with all defines
  input.clear();
  parser.parse_file(_inputfile, input);

  it = input.submodels_begin("Simulation");
  if (it == input.submodels_end("Simulation"))
    throw InitFailedException("\'Simulation\' block missing in input file.");
  const ModelOptions& global_opts = it->second;


  ModelOptions device_opts;
  device_opts["output_path"] = global_opts.get_option("resultpath", ".");
  it = input.submodels_begin("Device");
  device_opts += it->second;

  // create the device
  _device = Device::create(device_opts, dev_comm);

  // setup global options
  setup_globals(global_opts);


  // prepare the device
  _device->prepare();

  // create signal definitions
  it = input.submodels_begin("Signal");
  ModelOptions::submodel_iterator end = input.submodels_end("Signal");
  for ( ; it != end; ++it)
  {
    const ModelOptions& opts = it->second;
    SignalGenerator* sg = SignalGenerator::create(opts);
    sg->init();
  }


  // create and prepare modules
  it = input.submodels_begin("Module");
  end = input.submodels_end("Module");
  for ( ; it != end; ++it)
  {
    ModelOptions opts(it->second);

    // some global variables
    if (!opts.find_option("resultpath"))
      opts.set_option("resultpath", global_opts.get_option("resultpath", "."));

    if (!opts.find_option("scratchpath"))
      opts.set_option("scratchpath", global_opts.get_option("scratchpath",
          opts.get_option("resultpath", ".")));

    if (!opts.find_option("output_format"))
    {
      string default_format = "vtk";
      if (_device->get_mesh().mesh_dimension() == 1)
        default_format = "grace";
      opts.set_option("output_format",
          global_opts.get_option("output_format", default_format));
    }

    if (!opts.find_option("binary_output"))
      opts.set_option("binary_output",
          global_opts.get_option("binary_output", "true"));


    setup_module(_device, opts);
  }



  // initialize the device
  _device->init();




  // initialize the simulations, but only if they are not initialized yet
  auto simit = SimulationInterface::simulations_begin();
  const auto simend = SimulationInterface::simulations_end();
  for ( ; simit != simend; ++simit)
    if (!(*simit)->is_initialized())
      (*simit)->init();

}



void
Control::setup_globals(const ModelOptions& opts)
{
  using namespace boost::filesystem;

  _outputdir = opts.get_option("resultpath", _outputdir);

  // setup the logfile
  string logfile = opts.get_option("logfile", "");
  if (logfile.empty())
    logfile = _outputdir + "/" + Utils::basename(_inputfile) + ".log";



  // for the moment, if several processes try to use the same log file
  // we let everyone use a different one
  //if (TiberCad::get_mpi_comm().semiverify(&logfile))
  //{
  //  logfile = logfile + "." + to_string(TiberCad::get_mpi_comm().rank());
  //}
  logfile = logfile + "." + InputParser::get_defined("MPI_DEV_KEY");



  Messages::set_log_file(logfile, _device->get_communicator(), 0);

  Messages::info("Writing log to " + logfile);

  // Copy input file
  if (opts.get_option("backup_inputfile", true))
  {
    path in(_inputfile);
    path out(_outputdir + "/" + Utils::basename(_inputfile) + ".tib");
    if (!equivalent(in, out))
      copy_file(in, out, copy_options::overwrite_existing);
  }

  {
    ostringstream os;
    os << "tiberCAD release " << TiberCad::version_string(true)
       << " (" << TiberCad::arch_string() << ")\n"
       << TiberCad::compilation_system();
    Messages::info(os.str());

  }
  Messages::newline();
  Messages::info("Input file: " + _inputfile);
  Messages::newline();

  Database::set_search_path(opts.get_option("searchpath", ""));
  //opts.delete_option("searchpath");

  DLLoader::prepend_to_library_path(opts.get_option("modellibpath", "."));
  //opts.delete_option("modellibpath");

  _output_format = opts.get_option("output_format", "vtk");

  {
    Messages m;
    m.info("Initialize global simulation options");
    m.indent();

    // initialize global simulation options
    SimulationOptions::initialize(opts);

    std::string modellibpath;
    DLLoader::get_library_path(modellibpath);

    ostringstream os;
    os << "Simulation temperature: "
      << SimulationOptions::temperature << " K" << endl
      << "Database search path  : "
      << Database::get_search_path() << endl
      << "Module library path   : "
      << modellibpath << endl
      << "Log file              : "
      << logfile << endl;
    m.info(os.str());
    m.newline();
  }

  // create a global variable for time
  VariableValue::check_and_register("$time", _time);

  {
    int omp_procs = omp_get_max_threads();
    int mpi_comm = TiberCad::get_mpi_comm().size();

    if ((omp_procs > 1) || (mpi_comm) > 1)
    {
      Messages m;
      m.info("Parallelization");
      m.indent();

      ostringstream os;
      os << "MPI communicator size : "
         << mpi_comm << endl
         << "OMP threads           : "
         << omp_procs << endl;
      m.info(os.str());
      m.newline();

      if ((omp_procs > 1) && (mpi_comm) > 1)
      {
        Messages::info("This is a hybrid MPI+openMP run. "
            "Assure that OMP_NUM_THREADS is set correctly.");
        m.newline();
      }


      if (mpi_comm > 1)
      {
        char buffer[MPI_MAX_PROCESSOR_NAME];
        int buflen;
        MPI_Get_processor_name(buffer, &buflen);
        os.str("");
        os << "This is MPI rank " << TiberCad::get_mpi_comm().rank()
          << " out of " << mpi_comm
          << ", running on host " << buffer;
        Messages::info(os.str());
        Messages::newline();
      }

      m.newline();
    }
  }


  // what we want solve
  opts.get_option("solve", _solve_list);


  // create output directory
  path outpath(_outputdir);
  if (!exists(outpath))
  {
    // we catch any error here without doing anything yet
    try {
      create_directories(outpath);
    }
    catch (...) {}
  }

  if (!(exists(outpath) && is_directory(outpath)))
  {
    string msg("Cannot create or use '");
    msg += outpath.string() + "' as output directory.";
    throw InitFailedException(msg);
  }
}




void
Control::setup_module(Device* device, const ModelOptions& opts)
{
  string modelname = opts.get_name();

  IDSet phys_regions;
  string physreg = opts.get_option("regions", "all");
  device->extract_physical_regions(physreg, phys_regions);

  // get the user defined name (if defined...)
  string simulation_name = opts.get_option("name", modelname);

  Messages m;

  m.newline();
  m.info("Setting up simulation of type \'"
      + modelname + "\' ...");
  m.indent();


  //
  // now create the simulation
  //
  SimulationInterface* sim = SimulationInterface::create(modelname, opts);
  if (sim == NULL)
    throw ModelErrorException(
        "Unknown simulation type: " + modelname);


  sim->verbose() = SimulationOptions::verbose();

  sim->setup_environment(*device, phys_regions);

  sim->prepare();
}



double
Control::get_time(void)
{
  return _time;
}



void
Control::run_simulation(void)
{

  ostringstream os;
  os << "We solve: ";
  unsigned int n = _solve_list.size();
  for (unsigned int i = 0; i < n; i++)
    os << _solve_list[i] << " ";
  os << endl;
  Messages::info(os.str());

  vector<SimulationInterface*> simulations(n);

  // first check if we can find all simulations
  // We also let them solve the equilibrium
  for (unsigned int i = 0; i < n; i++)
  {
    SimulationInterface* sim = SimulationInterface::find_simulation(_solve_list[i]);

    if (sim == NULL)
      throw SolveFailedException("Simulation not found: " + _solve_list[i]);

    simulations[i] = sim;
  }

  Utils::Timer tt;

  // now run them
  for (unsigned int i = 0; i < n; i++)
  {
    SimulationInterface* sim = simulations[i];

    try
    {
      sim->solve();
      sim->plot();
    }
    catch (runtime_error& e)
    {
      // we plot anyway, maybe it helps to identify the problem
      //sim->plot();
      ostringstream s;
      s << "Solve of " << sim->get_name() << " failed." << endl <<
           "    Cause: " << e.what();
      throw SolveFailedException(s.str());
    }
    catch (...)
    {
      ostringstream s;
      s << "Solve of " << sim->get_name() << " failed for unknown reason.";
      throw SolveFailedException(s.str());
    }

  }
  os.str("");
  os << "Total solve time: " << tt.elapsed_string();
  Messages::newline();
  Messages::info(os.str());
}



void
Control::plot_all(void)
{
  SimulationInterface::SimulationIterator
    simit(SimulationInterface::simulations_begin());
  const SimulationInterface::SimulationIterator
    simend(SimulationInterface::simulations_end());

  for ( ; simit != simend; ++simit)
    (*simit)->plot();
}






