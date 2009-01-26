// $Id$


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
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationInterface.h"
#include "PetscRuntimeError.h"
#include "AtomisticStructure.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include <iostream>
#include <vector>
#include <set>


using namespace std;




Control::Control(const std::string& inputfile)
  : _inputfile(inputfile),
    _device(0),
    _database(0),
    _outputdir("."),
    _filename_suffix("")
{
  // we check here if the input file exists
  ifstream infile;
  infile.open(_inputfile.c_str());
  if (infile.fail() || !infile.good() || (infile.rdbuf()->in_avail() == 0))
  {
    infile.close();
    throw InitFailedException("Input file is invalid.");
  }

  infile.close();
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

  delete _database;
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

  _database = new Database();
  Material::set_database(*_database);

  // create the device, simulations and models
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
Control::create_device(void)
{
  using namespace boost::filesystem;

#ifdef DEBUG
  cerr << "Control::create_device() begin" << endl;
#endif

  InputParser parser(_inputfile);

  ModelOptions opts = parser.read_parameters("Simulation");

  _database->set_search_path(opts.get_option("searchpath", ""));
  opts.delete_option("searchpath");

  DLLoader::prepend_to_library_path(opts.get_option("modellibpath", "."));
  opts.delete_option("modellibpath");

#ifdef DEBUG
  cerr << " initialize global simulation options" << endl;
#endif
  // initialize global simulation options
  SimulationOptions::initialize(opts);


  //! read the variables we want to plot
  vector<string> vars;
  opts.get_option("plot", vars);
  opts.delete_option("plot");
  for (unsigned int i = 0; i < vars.size(); i++)
    _plotvariables.insert(vars[i]);

#ifdef DEBUG
  cerr << " create output directory" << endl;
#endif
  // create output directory
  _outputdir = opts.get_option("resultpath", _outputdir);
  opts.delete_option("resultpath");
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
    string msg("Cannot create ore use '");
    msg += outpath.string() + "' as output directory.";
    throw InitFailedException(msg);
  }


  _output_format = opts.get_option("output_format", "gmv");


#ifdef DEBUG
  cerr << " create device" << endl;
#endif
  _device = Device::create(opts);


  // tell the device who controls it
  _device->set_control(this);

  ostringstream os;
  os << "Simulation temperature = " <<
    SimulationOptions::temperature << " K" << endl;
  Messages::info(os.str());


#ifdef DEBUG
  cerr << "Control::create_device() end" << endl;
#endif
}



void
Control::create_materials(void)
{

#ifdef DEBUG
  cerr << "Control::create_materials() begin" << endl;
#endif

  assert(_device != NULL);

  InputParser parser(_inputfile);

  parser.read_device();

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

    const string& material = data.get_material_name();
    Material* mat = Material::create(material, data.get_options());
    _device->set_material(mat, region_ids, data.get_region_name());
  }

#ifdef DEBUG
  cerr << "Control::create_materials() end" << endl;
#endif
}



void
Control::create_atomistic_structures(void)
{

#ifdef DEBUG
  cerr << "Control::create_atomistic_structures() begin" << endl;
#endif

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

#ifdef DEBUG
  cerr << "Control::create_atomistic_structures() end" << endl;
#endif
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
      _device->get_region_ids(region_ids_str[i], tmp_id);
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

#ifdef DEBUG
  cerr << "Control::setup_models() begin" << endl;
#endif

  assert(_device != NULL);

  // a reference to the device
  Device& device = *_device;

  InputParser parser(_inputfile);

  // parse the models section
  const multimap<const string, ModelStructure*>& models = parser.read_models();

  // we loop over all simulations to setup the models
  multimap<const string, ModelStructure*>::const_iterator modit(models.begin());
  const multimap<const string,
        ModelStructure*>::const_iterator modend(models.end());

  for ( ; modit != modend; ++modit)
  {
    ModelStructure* model_str = modit->second;

    const ModelOptions& simopts = model_str->get_model_options();
    const string& modelname = model_str->get_model_name();

    //
    // extract the physical regions
    //

    set<ID> phys_regions;
    const string& physreg = simopts.get_option("physical_regions", "all");
    if (physreg == "all")
      phys_regions = _device->get_region_ids();
    else
    {
      // we have to get it as vector (for the moment at least)
      // we read them as strings as they could be region names
      vector<string> preg;
      simopts.get_option("physical_regions", preg);

      vector<ID> preg_ids;

      const set<ID>& regs = _device->get_region_ids();
      const set<ID>::const_iterator id_end(regs.end());
      unsigned int n = preg.size();
      for (unsigned int i = 0; i < n; i++)
      {
        // first check if it is a region name
        _device->get_region_ids(preg[i], preg_ids);
        if (preg_ids.size() != 0)
        {
          for (unsigned int j = 0; j < preg_ids.size(); j++)
            phys_regions.insert(preg_ids[j]);
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
          phys_regions.insert(id);
        }
      }
    }


    //
    // now create the simulation
    //

    // read options for this simulation (from Solver section)
    ModelOptions solveropts(parser.read_parameters("Solver", modelname));

    // get the user defined name (if defined...)
    const string& simulation_name = simopts.get_option("simulation_name", "");

    // read also section with user defined name as label
    if (!simulation_name.empty())
    {
      solveropts["name"] = simulation_name;

      if (simulation_name != modelname)
        solveropts += parser.read_parameters("Solver", simulation_name);
    }

    // we put the parameters in the $Solver section as submodel parameters
    // so we can hand them over in a cleaner way
    solveropts.add_submodel("$Solver", solveropts);


    // TODO we should change the handling of options, for now we just put the
    // two sections together
    solveropts += simopts;
    solveropts.delete_option("simulation_name");
    solveropts.delete_option("physical_regions");

    SimulationInterface* sim =
      SimulationInterface::create(modelname, solveropts);
    if (sim == NULL)
      throw ModelErrorException(
          "Unknown simulation type: " + modelname);

    sim->set_control(this);

    _simulations[sim->get_name()] = sim;

    // create the environment
    SimulationEnvironment* env = new SimulationEnvironment(device, phys_regions);
    _simulation_environments[sim] = env;
    sim->set_environment(env);

    sim->verbose() = SimulationOptions::verbose();


    //
    // now we have to create the models
    //

    // what main model should be used?
    ModelOptions physopts(parser.read_parameters("Physics", modelname));

    // read also section with user defined name as label
    if (!simulation_name.empty() && (simulation_name != modelname))
      physopts += parser.read_parameters("Physics", simulation_name);

    // parse the submodels
    // NOTE: different submodels could have the same identifier
    {
      multimap<const string, ModelOptions>& physmodels =
        model_str->get_physical_model_map();

      multimap<const string,
        ModelOptions>::iterator it(physmodels.begin());
      multimap<const string,
        ModelOptions>::iterator end(physmodels.end());
      for ( ; it != end; ++it)
      {
        // we set the name to the model type if not explicitly
        // given by user
        if (!(it->second).find_option("name"))
          (it->second)["name"] = it->first;

        physopts.add_submodel(it->first, it->second);
      }
    }


    // we have to do this for each material!
    set<ID>::iterator it(phys_regions.begin());
    const set<ID>::iterator end(phys_regions.end());
    for ( ; it != end; ++it)
    {
      Material* mat = device.get_material(*it);

      if (mat == NULL)
      {
        ostringstream s;
        s << "Physical region " << *it <<
          " has no material associated!";
        throw InitFailedException(s.str());
      }


      // here we actually create the model
      PhysicalModel* model = sim->create_physical_model(physopts, mat);

      // NOTE: model could be NULL, but we don't care about. Who tells us that
      // every simulation necessarily needs a model?
      if (model != NULL)
        mat->add_model(model, sim->get_id());
    }


    //
    // and now... the boundary conditions
    //

    map<ID, RegionStructure>& bc_map = model_str->get_model_BC_map();
    map<ID, RegionStructure>::iterator bdit(bc_map.begin());
    const map<ID, RegionStructure>::iterator bdend(bc_map.end());

    for ( ; bdit != bdend; ++bdit)
    {
      ID id = bdit->first;
      const RegionStructure& data = bdit->second;

      vector<ID> ids;
      Utils::extract_vector(data.get_region_ID(), ids);

      // if no numbers are specified we try to get them from the region name
      if (ids.size() == 0)
        _device->get_boundary_region_ids(data.get_region_name(), ids);
      else
        _device->set_boundary_region_name(data.get_region_name(), ids);

      if (ids.size() == 0)
      {
        ostringstream s;
        s << "Boundary region \'" << data.get_region_name() <<
          "\' is not consistent with mesh.";
        throw InitFailedException(s.str());
      }

      set<ID> region_ids;
      for (unsigned int i = 0; i < ids.size(); i++)
        region_ids.insert(ids[i]);

#ifdef DEBUG
      cerr << "Add boundary \'" << data.get_region_name()
        << "\' (region nr.";
      for (unsigned int i = 0; i < ids.size(); i++)
        cerr << " " << ids[i];
      cerr << ")\n";
#endif

      const ModelOptions& bdopts = data.get_options();

      Boundary* bd = new Boundary(data.get_region_name(), env, region_ids);
      bd->set_area_factor(bdopts.get_option("area_factor", 1.0));

      BoundaryProperties* bdprop = sim->create_boundary_model(bdopts);

      // NOTE: bdprop could be NULL, but we don't care about. Who tells us that
      // every simulation necessarily needs a boundary model?
      if (bdprop != NULL)
        bd->add_boundary_properties(bdprop, sim->get_id());

    }

    // prepare some of the environments internals (lists of elements etc.)
    env->prepare();

  } // end loop over simulations

  //
  // check for some special simulations (sweep, selfconsistency)
  //

  // first selfconsistency
  const multimap<string, ModelOptions>& scs =
    parser.read_subblocks("Solver", "Selfconsistent");

  multimap<string, ModelOptions>::const_iterator sc_it(scs.begin());
  multimap<string, ModelOptions>::const_iterator sc_end(scs.end());
  for ( ; sc_it != sc_end; ++sc_it)
  {
    const ModelOptions& solveropts = sc_it->second;
    if (!solveropts.is_empty())
    {
      SimulationInterface* sim =
        SimulationInterface::create("selfconsistent", solveropts);

      if (sim == NULL)
        throw ModelErrorException("Could not create Selfconsistent simulation");

      sim->set_control(this);
      sim->verbose() = 0;
      sim->set_name(sc_it->first);
      _simulations[sim->get_name()] = sim;
    }
  }

  if (scs.size() == 0)
  {
    const ModelOptions& solveropts =
      parser.read_parameters("Solver", "selfconsistent");
    if (!solveropts.is_empty())
    {
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

      Messages::warning("The definition of a selfconsistent simulation "
          "outside of a \'Selfconsistent\' block is deprecated.");
    }
  }



  // then the sweeps

  {
    // we might have several sweeps
    const multimap<string, ModelOptions>& sws =
      parser.read_subblocks("Solver", "Sweep");

    multimap<string, ModelOptions>::const_iterator sw_it(sws.begin());
    multimap<string, ModelOptions>::const_iterator sw_end(sws.end());
    for ( ; sw_it != sw_end; ++sw_it)
    {
      const ModelOptions& solveropts = sw_it->second;
      if (!solveropts.is_empty())
      {
        SimulationInterface* sim =
          SimulationInterface::create("sweep", solveropts);

        if (sim == NULL)
          throw ModelErrorException("Could not create Selfconsistent simulation");

        sim->set_control(this);
        sim->verbose() = 0;
        sim->set_name(sw_it->first);
        _simulations[sim->get_name()] = sim;
      }
    }

    // if no Sweep block is found, we look for other definitions to be compatible
    // with older TiberCAD version
    if (sws.size() == 0)
    {
      bool warning = false;
      ModelOptions sweepopts = parser.read_parameters("Solver", "sweep");
      if (!sweepopts.is_empty())
      {
        warning = true;
        if (!sweepopts.find_option("name"))
          sweepopts["name"] = "sweep";
        SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
        sim->set_control(this);
        //sim->verbose() = SimulationOptions::verbose();
        sim->verbose() = 0;
        _simulations[sim->get_name()] = sim;
      }

      sweepopts = parser.read_parameters("Solver", "sweep_1");
      if (!sweepopts.is_empty())
      {
        warning = true;
        if (!sweepopts.find_option("name"))
          sweepopts["name"] = "sweep_1";
        SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
        sim->set_control(this);
        //sim->verbose() = SimulationOptions::verbose();
        sim->verbose() = 0;
        _simulations[sim->get_name()] = sim;
      }

      sweepopts = parser.read_parameters("Solver", "sweep_2");
      if (!sweepopts.is_empty())
      {
        warning = true;
        if (!sweepopts.find_option("name"))
          sweepopts["name"] = "sweep_2";
        SimulationInterface* sim = SimulationInterface::create("sweep", sweepopts);
        sim->set_control(this);
        //sim->verbose() = SimulationOptions::verbose();
        sim->verbose() = 0;
        _simulations[sim->get_name()] = sim;
      }

      if (warning)
      {
        Messages::warning("The definition of sweeps outside of a \'Sweep\' "
            "block is deprecated.");
      }
    }
  }


#ifdef DEBUG
  cerr << "Control::setup_models() end" << endl;
#endif

}




  void
Control::run_simulation(void) throw (SolveFailedException)
{

  InputParser parser(_inputfile);

  const ModelOptions& opts = parser.read_parameters("Simulation");

  vector<string> names;
  opts.get_option("solve", names);

  ostringstream os;
  os << "We solve: ";
  unsigned int n = names.size();
  for (unsigned int i = 0; i < n; i++)
    os << names[i] << " ";
  os << endl;
  Messages::info(os.str());

  vector<SimulationInterface*> simulations(n);

  // first check if we can find all simulations
  // We also let them solve the equilibrium
  for (unsigned int i = 0; i < n; i++)
  {
    SimulationInterface* sim = find_simulation(names[i]);

    if (sim == NULL)
      throw SolveFailedException("Simulation not found: " + names[i]);

    simulations[i] = sim;
  }

  // now run them
  for (unsigned int i = 0; i < n; i++)
  {
    SimulationInterface* sim = simulations[i];

    // Let's do this in the models
    /*
    // we try first to solve the equilibrium solution, which
    // could be useful for most simulations
    try
    {
      sim->solve_equilibrium();
    }
    catch (runtime_error& e)
    {
      ostringstream s;
      s << "Control: Problem in solving equilibrium." << endl <<
        "    Cause: " << e.what();
      throw SolveFailedException(s.str());
    }
    catch (...)
    {
      ostringstream s;
      s << "Control: Solve of " << sim->get_name() << " failed." << endl <<
        "    Cause: Unknown";
      throw SolveFailedException(s.str());
    }
    */

    // now the actual solve
    try
    {
      sim->solve();
      sim->plot();
    }
    catch (runtime_error& e)
    {
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



