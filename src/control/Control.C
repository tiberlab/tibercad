// $Id$


#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "InputParser.h"
#include "Control.h"
#include "Database.h"
#include "DLLoader.h"
#include "Utils.h"
#include "Variable.h"
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
  : _inputfile(""),
    _device(0),
    _outputdir(".")
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

  Device::destroy(_device);

  // clear all variables
  VariableValue::clear_all();

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

  _device = Device::create(it->second);
  delete _device;

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
  _device = Device::create(device_opts);

  // setup global options
  setup_globals(global_opts);

  // prepare the device
  _device->prepare();


  // create and prepare modules
  it = input.submodels_begin("Module");
  ModelOptions::submodel_iterator end = input.submodels_end("Module");
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
      opts.set_option("output_format", global_opts.get_option("output_format", "vtk"));

    if (!opts.find_option("binary_output"))
      opts.set_option("binary_output", global_opts.get_option("binary_output", "true"));


    setup_module(_device, opts);
  }



  // initialize the device
  _device->init();




  // initialize the simulations, but only if they are not initialized yet
  SimulationInterface::SimulationIterator simit(SimulationInterface::simulations_begin());
  const SimulationInterface::SimulationIterator simend(SimulationInterface::simulations_end());
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


  Messages::info("Writing log to " + logfile);

  Messages::set_log_file(logfile, _device->get_communicator(), 0);

  {
    ostringstream os;
    os << "tiberCAD release " << TiberCad::version_string()
          << " ("
          << "rev. " << TiberCad::software_revision() << ", "
          << TiberCad::arch_string() << ")";
    Messages::info(os.str());

    if (TiberCad::get_mpi_comm().size() > 1)
    {
      char buffer[MPI_MAX_PROCESSOR_NAME];
      int buflen;
      MPI_Get_processor_name(buffer, &buflen);
      os.str("");
      os << "MPI: rank " << TiberCad::get_mpi_comm().rank()
            << " out of " << TiberCad::get_mpi_comm().size()
            << ", running on host " << buffer;
      Messages::newline();
      Messages::info(os.str());
    }
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

    ostringstream os;
    os << "Simulation temperature: "
      << SimulationOptions::temperature << " K" << endl
      << "Database search path  : "
      << Database::get_search_path() << endl
      //<< "Output directory      : "
      //<< _outputdir << endl
      << "Log file              : "
      << logfile << endl
      //<< "Output file format    : "
      //<< _output_format << endl
      << endl;
    m.info(os.str());
    m.newline();
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

  // the following stuff does not need to be done if the simulation is a task
  if (!sim->is_task())
  {

    SimulationEnvironment* env = NULL;
    sim->setup_environment(*device, phys_regions);
    env = &sim->get_environment();




    //
    // and now... the boundary conditions
    //

    // we accept the following keywords for boundaries:
    //   Contact, Boundary, Interface
    vector<string> keys(4);
    keys[0] = "Contact";
    keys[1] = "Boundary";
    keys[2] = "Interface";
    keys[3] = "BoundaryCondition";

    for (size_t i = 0; i < keys.size(); ++i)
    {
      ModelOptions::const_submodel_iterator it = opts.submodels_begin(keys[i]);
      const ModelOptions::const_submodel_iterator end = opts.submodels_end(keys[i]);
      for ( ; it != end; ++it)
        create_boundary(sim, it->second);

      // we remove them, so in the following we have only models
      //physopts.delete_submodels(keys[i]);
    }




    // the bulk physical models, all together
    ModelOptions bulk_opts;


    //
    // Next, we create all lower dimensional submodels
    // At the same time, we put all bulk Physics blocks together
    //
    // NOTE: only definitions with correct space dimensions will be
    //       added to the different PhysicalObject instances

    ModelOptions::const_submodel_iterator ph_it = opts.submodels_begin("Physics");
    const ModelOptions::const_submodel_iterator ph_end = opts.submodels_end("Physics");
    for ( ; ph_it != ph_end; ++ph_it)
    {
      ModelOptions physopts = ph_it->second;

      ModelOptions::submodel_iterator it = physopts.submodels_begin();
      const ModelOptions::submodel_iterator end = physopts.submodels_end();
      while (it != end)
      {
        ModelOptions::submodel_iterator tmpit(it);
        ++it;

        const ModelOptions& bdopts = tmpit->second;

        string physreg = bdopts.get_option("regions", "all");

        // if "all", it cannot be a lower dimensional region
        if (physreg == "all")
          continue;

        bool add = false;

        // It might be some lower dim model

        vector<string> ids_strings;
        Utils::extract_vector(physreg, ids_strings);

        for (unsigned int i = 0; i < ids_strings.size(); i++)
        {
          vector<ID> region_ids;
          device->get_boundary_region_ids(ids_strings[i], region_ids);

          // if it is no boundary region, we continue to the next region name
          if (region_ids.size() == 0)
          {
            // it could be the name of the Boundary object
            Boundary* bnd = env->get_boundary(ids_strings[i]);
            if (bnd != NULL)
              bnd->get_region_ids(region_ids);
          }

          if (region_ids.size() == 0)
            continue;

          // now it must be a lower dim model
          add = true;

          // NOTE: the model will only be added if a corresponding boundary ID
          // has been found
          for (unsigned int reg = 0; reg < region_ids.size(); reg++)
          {
            ID id = region_ids[reg];

            MaterialBoundary* bd;
            if ((bd = device->get_boundary_object(id)) != NULL)
            {
              PhysicalModel* pm = bd->get_model(sim->get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = sim->new_boundary_model(ModelOptions(), bd);
                bd->add_model(pm, sim->get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }

            EdgeObject* eo;
            if ((eo = device->get_edge_object(id)) != NULL)
            {
              PhysicalModel* pm = eo->get_model(sim->get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = sim->new_edge_model(ModelOptions(), eo);
                eo->add_model(pm, sim->get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }

            NodeObject* no;
            if ((no = device->get_node_object(id)) != NULL)
            {
              PhysicalModel* pm = no->get_model(sim->get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = sim->new_node_model(ModelOptions(), no);
                no->add_model(pm, sim->get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }
          }
        }

        // remove it from the map
        // NOTE: this means, we cannot specify the same model for bulk
        //       boundaries in a single block
        // TODO: maybe this can be relaxed?
        if (add)
          physopts.delete_submodel(tmpit);
      }

      // now add the remaining content to the common Physics options
      bulk_opts += physopts;
      bulk_opts.set_key((ph_it->second).get_key());
    }
    m.unindent();

    
    //
    // now we have to create the bulk models
    //
    //m.newline();
    m.info("Creating bulk models... ");
    m.indent();

    // we have to do this for each material!
    IDSet::iterator reg_it(phys_regions.begin());
    const IDSet::iterator reg_end(phys_regions.end());
    for ( ; reg_it != reg_end; ++reg_it)
    {
      ID reg_id = *reg_it;
      Material* mat = device->get_material(reg_id);

      m.info("Region: "+device->get_region_name(reg_id));
      if (mat == NULL)
      {
        ostringstream s;
        s << "Physical region " << device->get_region_name(reg_id) <<
            " has no material associated!";
        Messages::warning(s.str());
        continue;
        //throw InitFailedException(s.str());
      }

      // we only continue if the model has not already been added
      // this is important as a material can be assigned to different
      // regions
      // TODO (Is this true ??)
      if (mat->get_model(sim->get_id()) == NULL)
      {

        // the crystal structure
        //string crystal_structure(mat->get_structure());

        // that's not elegant, but it works
        ModelOptions opts(bulk_opts);
        opts.delete_all_submodels();

        // we add the crystal structure for bulk materials as this could
        // lead to different model implementations
        //opts["crystal_structure"] = crystal_structure;

        //
        // we parse the submodels for each region as they could be associated
        // one-by-one
        // NOTE: different submodels could have the same identifier
        //
        {
          ModelOptions::submodel_iterator it = bulk_opts.submodels_begin();
          const ModelOptions::submodel_iterator end = bulk_opts.submodels_end();
          for ( ; it != end; ++it)
          {
            bool add = true;

            ModelOptions modopts(it->second);

            // we have to check if it should be built for the current region
            // TODO to not allow for errors

            IDSet regs;
            string physreg = modopts.get_option("regions", "all");
            device->extract_physical_regions(physreg, regs);

            if (regs.count(reg_id) == 0) add = false;


            if (add)
            {
              // we add the crystal structure for bulk materials as this could
              // lead to different model implementations
              //modopts.set_option("crystal_structure", crystal_structure);

              // we set the name to the model type if not explicitly
              // given by user
              //if (!(mapit->second).find_option("name"))
              //  (mapit->second)["name"] = mapit->first;
              opts.add_submodel(it->first, modopts);
            }
          }
        }


        // here we actually create the model
        PhysicalModel* model = sim->new_bulk_model(opts, mat);

        // NOTE: model could be NULL, but we don't care about. Who tells us that
        // every simulation necessarily needs a model?
        mat->add_model(model, sim->get_id());
      }
    }
    m.unindent();

  }

  sim->prepare();
}





void
Control::create_boundary(SimulationInterface* sim, const ModelOptions& opts)
{
  if (sim->is_task()) return;

  SimulationEnvironment* env = &(sim->get_environment());

  if (env == NULL) return;

  Device& device = env->get_device();

  // first get region names
  string id_str = opts.get_option("regions", "");
  vector<string> ids_strings;
  Utils::extract_vector(id_str, ids_strings);

  string boundary_name = opts.get_name();

  // for the numeric IDs
  vector<ID> ids;

  unsigned int n_ids = ids_strings.size();
  // if no numbers are specified we try to get them from the region name
  if (n_ids == 0)
    device.get_boundary_region_ids(boundary_name, ids);
  else
  {
    vector<ID> tmp_id;
    for (unsigned int i = 0; i < n_ids; i++)
    {
      // either it is a name or a number
      // try first name
      device.get_boundary_region_ids(ids_strings[i], tmp_id);
      if (tmp_id.size() == 0)
      {
        ostringstream s;
        s << "Physical region \'" << ids_strings[i]
          << "\' (in boundary  \'" << boundary_name
          << "\') does not exist in mesh.";
        throw InitFailedException(s.str());
      }
      ids.insert(ids.end(), tmp_id.begin(), tmp_id.end());
    }
  }

  if (ids.size() == 0)
  {
    ostringstream s;
    s << "Boundary region \'" << boundary_name <<
        "\' is not consistent with mesh.";
    throw InitFailedException(s.str());
  }

  {
    ostringstream os;
    os << "Adding boundary \'" << boundary_name << "\'";
    Messages::info(os.str());
  }


  Boundary* bnd = env->get_boundary(boundary_name);
  if (bnd == NULL)
  {
    bnd = new Boundary(boundary_name, opts);
    bnd->set_region_ids(ids);
    env->add_boundary(bnd);

    //
    // this is the old way -->

    //BoundaryProperties* bdprop = sim->new_boundary_model(opts);

    // NOTE: bdprop could be NULL, but we don't care about. Who tells us that
    // every simulation necessarily needs a boundary model?
    //if (bdprop != NULL)
    //  bnd->add_boundary_properties(bdprop, sim->get_id());

    // <-- end of old way
    //
  }


  for (unsigned int i = 0; i < ids.size(); i++)
  {

    bool found = false;

    MaterialBoundary* bd;
    if ((bd = device.get_boundary_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = bd->get_model(sim->get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << sim->get_name();
        throw InitFailedException(os.str());
      }

      ID ida = bd->get_id_A();
      ID idb = bd->get_id_B();
      if (!sim->includes_region(ida) && !sim->includes_region(idb))
      {
        continue;

        Messages::warning("At least part of boundary \'" + boundary_name +
            "\' does not touch any region of simulation \'" +
            sim->get_name() + "\'");
      }

      pm = sim->new_boundary_model(opts, bd);
      bd->add_model(pm, sim->get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    EdgeObject* eo;
    if ((eo = device.get_edge_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = eo->get_model(sim->get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << sim->get_name();
        throw InitFailedException(os.str());
      }

      pm = sim->new_edge_model(opts, eo);
      eo->add_model(pm, sim->get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    NodeObject* no;
    if ((no = device.get_node_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = no->get_model(sim->get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << sim->get_name();
        throw InitFailedException(os.str());
      }

      pm = sim->new_node_model(opts, no);
      no->add_model(pm, sim->get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    if (!found)
    {
      ostringstream os;
      os << "Boundary \'" << boundary_name
                  << "\' does not exist.";
      throw InitFailedException(os.str());
    }

  }

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






