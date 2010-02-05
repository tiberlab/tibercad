// $Id$


#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "InputParser.h"
#include "RegionStructure.h"
#include "Control.h"
#include "Database.h"
#include "DLLoader.h"
#include "Utils.h"
#include "Messages.h"
#include "SimulationOptions.h"
#include "Device.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationInterface.h"
#include "PetscRuntimeError.h"
#include "AtomisticStructure.h"

#include <sstream>
#include <vector>
#include <set>

#include <csignal>


using namespace std;




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



Control::Control(void)
  : _inputfile(""),
    _device(0),
    _outputdir(".")
{

}



Control::~Control(void)
{
  simulation_iterator simit(_simulations.begin());
  const simulation_iterator simend(_simulations.end());

  for ( ; simit != simend; ++simit)
    SimulationInterface::destroy(*simit);

  _simulations.clear();

  EnvironmentMap::iterator envit(_simulation_environments.begin());
  const EnvironmentMap::iterator envend(_simulation_environments.end());
  for ( ; envit != envend; ++envit)
    delete envit->second;

  _simulation_environments.clear();

  Device::destroy(_device);

}


void
Control::invalidate_environments(void)
{
  EnvironmentMap::iterator envit(_simulation_environments.begin());
  const EnvironmentMap::iterator envend(_simulation_environments.end());
  for ( ; envit != envend; ++envit)
    envit->second->invalidate();
}



void
Control::init(void) throw (InitFailedException,
    ModelErrorException, DatabaseException)
{

  // TODO this is currently broken
  // setup the signal handler
  //SignalHandler::set_control(this);

  // we want to intercept SIGINT (Ctrl-C)
  //SignalHandler::activate_sigint();


  // we check here if the input file exists
  ifstream infile;
  infile.open(_inputfile.c_str());
  if (infile.fail() || !infile.good())
  {
    infile.close();
    throw InitFailedException("Input file is invalid.");
  }
  infile.close();


  // create the device, simulations and models
  setup_globals();
  create_device();
  create_materials();
  setup_clusters();
  create_atomistic_structures();
  setup_models();

  // initialize the device
  _device->init();


  // initialize the simulation environments
  EnvironmentMap::iterator envit(_simulation_environments.begin());
  const EnvironmentMap::iterator envend(_simulation_environments.end());
  for ( ; envit != envend; ++envit)
    envit->second->init();


  // initialize the simulations, but only if they are not initialized yet
  // (the latter should not happen, however)
  simulation_iterator simit(_simulations.begin());
  const simulation_iterator simend(_simulations.end());
  for ( ; simit != simend; ++simit)
    if (!(*simit)->is_initialized())
      (*simit)->init();

}



void
Control::setup_globals(void)
{
  using namespace boost::filesystem;

  InputParser parser(_inputfile);

  ModelOptions opts;
  parser.get_simulation_options(opts);

  // setup the logfile
  string logfile(Utils::basename(_inputfile) + ".log");
  logfile = opts.get_option("logfile", logfile);
  opts.delete_option("logfile");
  Messages::set_log_file(logfile);

  Messages::info("Input file: " + _inputfile);
  Messages::newline();

  Database::set_search_path(opts.get_option("searchpath", ""));
  opts.delete_option("searchpath");

  DLLoader::prepend_to_library_path(opts.get_option("modellibpath", "."));
  opts.delete_option("modellibpath");


  _outputdir = opts.get_option("resultpath", _outputdir);
  opts.delete_option("resultpath");

  _output_format = opts.get_option("output_format", "gmv");

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
      << "Output directory      : "
      << _outputdir << endl
      << "Log file              : "
      << logfile << endl
      << "Output file format    : "
      << _output_format << endl << endl;
    m.info(os.str());
    m.newline();
  }


  // what we want solve
  opts.get_option("solve", _solve_list);

  // read the variables we want to plot
  vector<string> vars;
  opts.get_option("plot", vars);
  opts.delete_option("plot");
  for (unsigned int i = 0; i < vars.size(); i++)
    _plotvariables.insert(vars[i]);


  // create output directory
  path outpath(_outputdir, native);
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
Control::create_device(void)
{
  Messages::debug("Control::create_device() begin");

  InputParser parser(_inputfile);

  ModelOptions opts;
  parser.get_simulation_options(opts);
  if (opts.find_option("mesh_units"))
  {
    Messages::warning("\"mesh_units\" in $Simulation section is deprecated "
        " and should be moved to the device options.");
  }
  if (opts.find_option("dimension"))
  {
    Messages::warning("\"dimension\" in $Simulation section is deprecated "
        " and should be moved to the device options.");
  }
  if (opts.find_option("meshfile"))
  {
    Messages::warning("\"meshfile\" in $Simulation section is deprecated "
        " and should be moved to the device options.");
  }
  if (opts.find_option("symmetry"))
  {
    Messages::warning("\"symmetry\" in $Simulation section is deprecated "
        " and should be moved to the device options.");
  }

  parser.read_device();
  opts += parser.get_device_options();

  // we pass the remaining options to the device
  _device = Device::create(opts);


  // tell the device who controls it
  _device->set_control(this);

  Messages::debug("Control::create_device() end");
}



void
Control::create_materials(void)
{

  Messages::debug("Control::create_materials() begin");

  assert(_device != NULL);

  InputParser parser(_inputfile);

  parser.read_device();

  Messages m;
  m.info("Create materials ...");
  m.indent();

  //
  // first we process the physical regions
  //
  const map<ID, RegionStructure>& device_map = parser.get_device_map();

  // iterate the map and create the materials
  map<ID, RegionStructure>::const_iterator mapit(device_map.begin());
  const map<ID, RegionStructure>::const_iterator mapend(device_map.end());
  for ( ; mapit != mapend; ++mapit)
  {
    const RegionStructure& data = mapit->second;

    // we read the region numbers as strings as they could be region names
    vector<string> region_ids_str;
    Utils::extract_vector(data.get_region_ID(), region_ids_str);

    // for the numeric region IDs
    vector<ID> region_ids;

    unsigned int n_ids = region_ids_str.size();
    // if no numbers are specified we try to get them from the region name
    if (n_ids == 0)
      _device->get_mesh_region_ids(data.get_region_name(), region_ids);
    else
    {
      vector<ID> tmp_id;
      for (unsigned int i = 0; i < n_ids; i++)
      {
        // either it is a name or a number
        // try first name
        _device->get_mesh_region_ids(region_ids_str[i], tmp_id);
        if (tmp_id.size() > 0)
          region_ids.insert(region_ids.end(), tmp_id.begin(), tmp_id.end());
        else
          region_ids.push_back(Utils::convert<unsigned int>(region_ids_str[i]));
      }
    }


    if (region_ids.size() == 0)
    {
      ostringstream s;
      s << "Physical region \'" << data.get_region_name() <<
        "\' is not consistent with mesh.";
      throw InitFailedException(s.str());
    }

    // get the common device options
    ModelOptions opts(parser.get_device_options());
    // The default material is Si
    // backwards compatibility
    string material = opts.get_option("mat", "Si");
    material = opts.get_option("material", material);

    opts += data.get_options();

    // backwards compatibility
    material = data.get_options().get_option("mat", material);
    material = data.get_options().get_option("material", material);
    opts.delete_option("mat");
    opts["material"] = material;

    Material* mat = Material::create(material, opts);
    _device->set_material(mat, region_ids, data.get_region_name());
  }

  m.unindent();
  m.info("Creation of materials done.");

  Messages::debug("Control::create_materials() end");
}



void
Control::create_atomistic_structures(void)
{

  Messages::debug("Control::create_atomistic_structures() begin");

  assert(_device != NULL);

  InputParser parser(_inputfile);

  parser.read_scale();

  //
  // and now we look for atomistic structures
  //
  const map<ID, RegionStructure>& atomistic_map = parser.get_atomistic_map();

  ModelOptions atomistic_options;

 // iterate the map and create the structures
  map<ID, RegionStructure>::const_iterator mapit(atomistic_map.begin());
  const map<ID, RegionStructure>::const_iterator mapend(atomistic_map.end());
  for ( ; mapit != mapend; ++mapit)
  {
    const RegionStructure& data = mapit->second;

    ModelOptions atomistic_options = data.get_options();

    const string& st_name = data.get_region_name();
    AtomisticStructure* st = AtomisticStructure::create(st_name, data.get_options());

    //Set a reference to device. It's needed for keeping track of region informations
    st->set_device( _device );

    //WARNING: For debugging purposes, initialization of
    //atomistic structures is here, but it's not the right place! (maybe it is...)
    st->init();

    // Defined atomistic structure is put in the atomistic_structure_map
    _device->set_atomistic_structure(st_name, st);

  }

  Messages::debug("Control::create_atomistic_structures() end");
}



void
Control::setup_clusters(void)
{
  InputParser parser(_inputfile);

  parser.read_device();

  const map<ID, RegionStructure>& cluster_map = parser.get_cluster_map();

  map<ID, RegionStructure>::const_iterator mapit(cluster_map.begin());
  const map<ID, RegionStructure>::const_iterator mapend(cluster_map.end());
  for ( ; mapit != mapend; ++mapit)
  {
    const RegionStructure& data = mapit->second;

    // we read the region numbers as strings as they could be region names
    vector<string> region_ids_str;
    Utils::extract_vector(data.get_region_ID(), region_ids_str);

    // for the numeric region IDs
    vector<ID> region_ids;

    unsigned int n_ids = region_ids_str.size();
    vector<ID> tmp_id;
    for (unsigned int i = 0; i < n_ids; i++)
    {
      // either it is a name or a number
      // try first name
      _device->get_active_region_ids(region_ids_str[i], tmp_id);
      if (tmp_id.size() > 0)
        region_ids.insert(region_ids.end(), tmp_id.begin(), tmp_id.end());
      else
        region_ids.push_back(Utils::convert<unsigned int>(region_ids_str[i]));
    }

    if (region_ids.size() > 0)
    {
      ostringstream os;
      os << "Setting up Cluster \'" << data.get_region_name()
        << "\' containing regions " << region_ids[0];
      for (size_t i = 1; i < region_ids.size(); i++)
        os << ", " << region_ids[i];
      Messages::info(os.str());

      _device->set_cluster(data.get_region_name(), region_ids);
    }
    else
      Messages::warning("Cluster \'" + data.get_region_name() + "\' is empty.");
  }
}





void
Control::setup_models(void) throw (InitFailedException, ModelErrorException)
{

  Messages::debug("Control::setup_models() begin");

  assert(_device != NULL);

  // a reference to the device
  Device& device = *_device;

  InputParser parser(_inputfile);

  typedef multimap<const string, ModelStructure*> ModelsMap;
  typedef map<string, ModelOptions> OptionsMap;

  // parse the models section
  const ModelsMap& models = parser.read_models();

  // get the blocks of the Physics section
  OptionsMap physics_opts;
  parser.get_physics_options_map(physics_opts);

  // get the blocks of the Solver section
  OptionsMap solver_opts;
  parser.get_solver_options_map(solver_opts);


  // we loop over all simulations to setup the models
  ModelsMap::const_iterator modit(models.begin());
  const ModelsMap::const_iterator modend(models.end());

  for ( ; modit != modend; ++modit)
  {
    ModelStructure* model_str = modit->second;

    ModelOptions simopts(model_str->get_model_options());
    const string& modelname = model_str->get_model_name();

    //
    // extract the physical regions
    //

    IDSet phys_regions;
    const string& physreg = simopts.get_option("physical_regions", "all");
    extract_physical_regions(physreg, phys_regions);

    // get the user defined name (if defined...)
    string simulation_name = simopts.get_option("simulation_name", "");
    simulation_name = simopts.get_option("name", simulation_name);

    // some cleanup
    if (!simulation_name.empty())
      simopts["name"] = simulation_name;
    simopts.delete_option("simulation_name");
    simopts.delete_option("physical_regions");

    Messages m;

    //
    // now create the simulation
    //
    m.newline();
    m.info("Setting up simulation of type \'"
        + modelname + "\' ...");
    m.indent();

    // read solver options for this simulation (from Solver section)
    ModelOptions solveropts;

    OptionsMap::iterator map_it(solver_opts.find(modelname));
    if (map_it != solver_opts.end())
      solveropts += map_it->second;

    // read also section with user defined name as label
    if (!simulation_name.empty())
    {
      map_it = solver_opts.find(simulation_name);
      if ((simulation_name != modelname) && (map_it != solver_opts.end()))
        solveropts += map_it->second;
    }

    // the main physical model -> PhysicalModel
    // we need this below
    ModelOptions physopts(simopts);


    // we put the parameters in the $Solver section as submodel parameters
    // so we can hand them over in a cleaner way
    simopts.add_submodel("$Solver", solveropts);
    // for compatibility we add also the solver options
    simopts += solveropts;


    SimulationInterface* sim =
      SimulationInterface::create(modelname, simopts);
    if (sim == NULL)
      throw ModelErrorException(
          "Unknown simulation type: " + modelname);

    sim->set_control(this);

    _simulations[sim->get_name()] = sim;

    // create the environment
    SimulationEnvironment* env =
      new SimulationEnvironment(device, phys_regions);
    _simulation_environments[sim] = env;
    sim->set_environment(env);

    sim->verbose() = SimulationOptions::verbose();


    //
    // now we have to create the models
    //
    m.info("Creating physical models... ", false);

    // TODO only for backwards compatibility
    map_it = physics_opts.find(modelname);
    if (map_it != physics_opts.end())
    {
      physopts += map_it->second;
      Messages::warning("The $Physics section is deprecated. \nOptions should be put"
          " into the \'options\' block of the model instead.");
    }

    // read also section with user defined name as label
    map_it = physics_opts.find(simulation_name);
    if (!simulation_name.empty() && (simulation_name != modelname) &&
        (map_it != physics_opts.end()))
    {
      physopts += map_it->second;
      Messages::warning("The $Physics section is deprecated. Options should be put"
          " into the \'options\' block of the model instead.");
    }


    // we have to do this for each material!
    IDSet::iterator it(phys_regions.begin());
    const IDSet::iterator end(phys_regions.end());
    for ( ; it != end; ++it)
    {
      ID reg_id = *it;
      Material* mat = device.get_material(reg_id);

      if (mat == NULL)
      {
        ostringstream s;
        s << "Physical region " << reg_id <<
          " has no material associated!";
        Messages::warning(s.str());
        continue;
        //throw InitFailedException(s.str());
      }

      // we only continue if the model has not already been added
      // this is important as a material can be assigned to different
      // regions
      if (mat->get_model(sim->get_id()) == NULL)
      {

        // the crystal structure
        string crystal_structure(mat->get_structure());

        ModelOptions opts(physopts);
        // we add the crystal structure for bulk materials as this could
        // lead to different model implementations
        opts["crystal_structure"] = crystal_structure;

        //
        // we parse the submodels for each region as they could be associated
        // one-by-one
        // NOTE: different submodels could have the same identifier
        //
        {
          multimap<const string, ModelOptions>& physmodels =
            model_str->get_physical_model_map();

          multimap<const string,
            ModelOptions>::iterator mapit(physmodels.begin());
          multimap<const string,
            ModelOptions>::iterator mapend(physmodels.end());
          for ( ; mapit != mapend; ++mapit)
          {
            bool add = true;
            // we have to check if it should be built for the current region
            if ((mapit->second).find_option("restrict_to_region"))
            {
              IDSet regs;
              const string& physreg =
                (mapit->second).get_option("restrict_to_region", "all");
              extract_physical_regions(physreg, regs);
              if (regs.count(reg_id) == 0) add = false;
            }


            if (add)
            {
              // we add the crystal structure for bulk materials as this could
              // lead to different model implementations
              (mapit->second)["crystal_structure"] = crystal_structure;

              // we set the name to the model type if not explicitly
              // given by user
              if (!(mapit->second).find_option("name"))
                (mapit->second)["name"] = mapit->first;
              opts.add_submodel(mapit->first, mapit->second);
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

    m.info("done");


    //
    // and now... the boundary conditions
    //
    m.info("Setup of boundary models...", false);

    map<ID, RegionStructure>& bc_map = model_str->get_model_BC_map();
    map<ID, RegionStructure>::iterator bdit(bc_map.begin());
    const map<ID, RegionStructure>::iterator bdend(bc_map.end());

    for ( ; bdit != bdend; ++bdit)
    {
      //ID id = bdit->first;
      const RegionStructure& data = bdit->second;

      // first get region names
      vector<string> ids_strings;
      Utils::extract_vector(data.get_region_ID(), ids_strings);

      // for the numeric IDs
      vector<ID> ids;

      unsigned int n_ids = ids_strings.size();
      // if no numbers are specified we try to get them from the region name
      if (n_ids == 0)
        _device->get_boundary_region_ids(data.get_region_name(), ids);
      else
      {
        vector<ID> tmp_id;
        for (unsigned int i = 0; i < n_ids; i++)
        {
          // either it is a name or a number
          // try first name
          _device->get_boundary_region_ids(ids_strings[i], tmp_id);
          if (tmp_id.size() > 0)
            ids.insert(ids.end(), tmp_id.begin(), tmp_id.end());
          else
            ids.push_back(Utils::convert<unsigned int>(ids_strings[i]));
        }
      }

      if (ids.size() == 0)
      {
        ostringstream s;
        s << "Boundary region \'" << data.get_region_name() <<
          "\' is not consistent with mesh.";
        throw InitFailedException(s.str());
      }

      IDSet region_ids;
      for (unsigned int i = 0; i < ids.size(); i++)
        region_ids.insert(ids[i]);

      {
        ostringstream os;
        os << "Add boundary \'" << data.get_region_name()
          << "\' (region nr.";
        for (unsigned int i = 0; i < ids.size(); i++)
          os << " " << ids[i];
        os << ")";
        Messages::debug(os.str());
      }

      const ModelOptions& bdopts = data.get_options();

      //
      // this is the old way -->

      Boundary* bd = new Boundary(data.get_region_name(), env, region_ids);
      bd->set_area_factor(bdopts.get_option("area_factor", 1.0));

      BoundaryProperties* bdprop = sim->new_boundary_model(bdopts);

      // NOTE: bdprop could be NULL, but we don't care about. Who tells us that
      // every simulation necessarily needs a boundary model?
      if (bdprop != NULL)
        bd->add_boundary_properties(bdprop, sim->get_id());

      // <-- end of old way
      //

      for (IDSet::const_iterator it(region_ids.begin());
           it != region_ids.end(); ++it)
      {
        ID id = *it;

        MaterialBoundary* bd;
        if ((bd = device.get_boundary_object(id)) != NULL)
        {
          Material* matA = bd->get_material_A();
          Material* matB = bd->get_material_B();
          PhysicalModel* pm = sim->new_boundary_model(bdopts, matA, matB);
          bd->add_model(pm, sim->get_id());
        }

        EdgeObject* eo;
        if ((eo = device.get_edge_object(id)) != NULL)
        {
          PhysicalModel* pm = sim->new_edge_model(bdopts);
          eo->add_model(pm, sim->get_id());
        }

        NodeObject* no;
        if ((no = device.get_node_object(id)) != NULL)
        {
          PhysicalModel* pm = sim->new_node_model(bdopts);
          no->add_model(pm, sim->get_id());
        }
      }

    }
    m.info("done");

    // prepare some of the environments internals (lists of elements etc.)
    env->prepare();

  } // end loop over simulations

  //
  // check for some special simulations (sweep, selfconsistency)
  //

  // first selfconsistency
  OptionsMap::iterator map_it(solver_opts.find("Selfconsistent"));
  if (map_it != solver_opts.end())
  {
    ModelOptions& sc_opts = map_it->second;

    ModelOptions::const_submodel_iterator sc_it(sc_opts.submodels_begin());
    const ModelOptions::const_submodel_iterator sc_end(sc_opts.submodels_end());

    for ( ; sc_it != sc_end; ++sc_it)
    {
      const ModelOptions& solveropts = sc_it->second;
      if (!solveropts.is_empty())
      {
        Messages m;
        m.newline();
        m.info("Setting up a selfconsistent simulation ("
            + sc_it->first + ") ...");
        m.indent();

        SimulationInterface* sim =
          SimulationInterface::create("selfconsistent", solveropts);

        if (sim == NULL)
          throw ModelErrorException("Could not create Selfconsistent simulation");

        sim->set_control(this);
        sim->verbose() = 0;
        sim->set_name(sc_it->first);
        _simulations[sim->get_name()] = sim;
        m.unindent();
      }
    }
  }

  map_it = solver_opts.find("selfconsistent");
  if (map_it != solver_opts.end())
  {
    const ModelOptions& solveropts = map_it->second;

    Messages::warning("The definition of a selfconsistent simulation "
        "outside of a \'Selfconsistent\' block is deprecated.");

    Messages m;
    m.newline();
    m.info("Setting up a selfconsistent simulation ("
        + map_it->first + ") ...");
    m.indent();

    SimulationInterface* sim =
      SimulationInterface::create("selfconsistent", solveropts);

    if (sim == NULL)
    {
      string msg("No such simulation type: selfconsistent (flavour: ");
      msg += solveropts.get_option("flavour", "");
      throw ModelErrorException(msg);
    }
    sim->set_control(this);
    //sim->verbose() = SimulationOptions::verbose();
    sim->verbose() = 0;
    _simulations[sim->get_name()] = sim;

    m.unindent();
  }



  // then the sweeps

  {
    OptionsMap::iterator map_it(solver_opts.find("Sweep"));
    if (map_it != solver_opts.end())
    {
      ModelOptions& sw_opts = map_it->second;

      ModelOptions::const_submodel_iterator sw_it(sw_opts.submodels_begin());
      const ModelOptions::const_submodel_iterator sw_end(sw_opts.submodels_end());


      for ( ; sw_it != sw_end; ++sw_it)
      {
        const ModelOptions& solveropts = sw_it->second;
        if (!solveropts.is_empty())
        {
          Messages m;
          m.newline();
          m.info("Setting up a parameter sweep ("
              + sw_it->first + ") ...");
          m.indent();

          SimulationInterface* sim =
            SimulationInterface::create("sweep", solveropts);

          if (sim == NULL)
            throw ModelErrorException("Could not create sweep simulation");

          sim->set_control(this);
          sim->verbose() = 0;
          sim->set_name(sw_it->first);
          _simulations[sim->get_name()] = sim;
          m.unindent();
        }
      }
    }

    // we look for other definitions to
    // be compatible with older TiberCAD version
    bool warning = false;

    // the following is allowed for ease of use
    map_it = solver_opts.find("sweep");
    if (map_it != solver_opts.end())
    {
      Messages m;
      m.newline();
      m.info("Setting up a parameter sweep ("
          + map_it->first + ") ...");
      m.indent();

      ModelOptions sweepopts(map_it->second);

      if (!sweepopts.find_option("name"))
        sweepopts["name"] = "sweep";
      SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
      sim->set_control(this);
      //sim->verbose() = SimulationOptions::verbose();
      sim->verbose() = 0;
      _simulations[sim->get_name()] = sim;
      m.unindent();
    }

    map_it = solver_opts.find("sweep_1");
    if (map_it != solver_opts.end())
    {
      Messages m;
      m.newline();
      m.info("Setting up a parameter sweep ("
          + map_it->first + ") ...");
      m.indent();

      ModelOptions sweepopts(map_it->second);

      warning = true;
      if (!sweepopts.find_option("name"))
        sweepopts["name"] = "sweep_1";
      SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
      sim->set_control(this);
      //sim->verbose() = SimulationOptions::verbose();
      sim->verbose() = 0;
      _simulations[sim->get_name()] = sim;
      m.unindent();
    }

    map_it = solver_opts.find("sweep_2");
    if (map_it != solver_opts.end())
    {
      Messages m;
      m.newline();
      m.info("Setting up a parameter sweep ("
          + map_it->first + ") ...");
      m.indent();

      ModelOptions sweepopts(map_it->second);

      warning = true;
      if (!sweepopts.find_option("name"))
        sweepopts["name"] = "sweep_2";
      SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
      sim->set_control(this);
      //sim->verbose() = SimulationOptions::verbose();
      sim->verbose() = 0;
      _simulations[sim->get_name()] = sim;
      m.unindent();
    }

    if (warning)
    {
      Messages::warning("The definition of sweeps outside of a \'Sweep\' "
          "block is deprecated.");
    }
  }

  Messages::debug("Control::setup_models() end");
}




void
Control::extract_physical_regions(const std::string& str, IDSet& ids)
{

  if (str == "all")
    ids = _device->get_active_region_ids();
  else
  {
    // we have to get it as vector (for the moment at least)
    // we read them as strings as they could be region names
    vector<string> preg;
    Utils::extract_vector(str, preg);

    vector<ID> preg_ids;

    const IDSet& regs = _device->get_active_region_ids();
    const IDSet::const_iterator id_end(regs.end());
    unsigned int n = preg.size();
    for (unsigned int i = 0; i < n; i++)
    {
      // first check if it is a region name
      _device->get_active_region_ids(preg[i], preg_ids);
      if (preg_ids.size() != 0)
      {
        for (unsigned int j = 0; j < preg_ids.size(); j++)
          ids.insert(preg_ids[j]);
      }
      else
      {
        // it has to be a region number
        ID id = Utils::convert<ID>(preg[i]);

        if (regs.find(id) == id_end)
        {
          ostringstream s;
          s << "Physical region " << id <<
            " does not exist in mesh file.";
          throw InitFailedException(s.str());
        }
        ids.insert(id);
      }
    }
  }
}




void
Control::run_simulation(void) throw (SolveFailedException)
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
    SimulationInterface* sim = find_simulation(_solve_list[i]);

    if (sim == NULL)
      throw SolveFailedException("Simulation not found: " + _solve_list[i]);

    simulations[i] = sim;
  }

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
      sim->plot();
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
}



void
Control::plot_all(void)
{
  simulation_iterator simit(_simulations.begin());
  const simulation_iterator simend(_simulations.end());

  for ( ; simit != simend; ++simit)
    (*simit)->plot();
}



SimulationInterface*
Control::find_simulation(const string& name) const
{
  SimulationInterface* sim = SimulationInterface::find_simulation(name);

  if (sim != NULL)
  {
    SimulationMap::const_iterator it = _simulations.find(sim->get_name());
    if (it == _simulations.end())
      sim = NULL;
  }

  return sim;
}



const std::string&
Control::get_filename_suffix(void) const
{
  _filename_suffix_str = "";
  list<string>::const_iterator it(_filename_suffix.begin());
  const list<string>::const_iterator end(_filename_suffix.end());
  for ( ; it != end; ++it)
    _filename_suffix_str += "_" + *it;

  return _filename_suffix_str;
}


