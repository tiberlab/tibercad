// $Id$


#include "InputParser.h"
#include "Control.h"
#include "Database.h"
#include "Utils.h"
#include "SimulationOptions.h"
#include "Device.h"
#include "Material.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationInterface.h"


#include <iostream>
#include <vector>
#include <set>


using namespace std;


Control::~Control(void)
{
  SimulationMap::iterator simit(_simulations.begin());
  const SimulationMap::iterator simend(_simulations.end());

  for ( ; simit != simend; ++simit)
    SimulationInterface::destroy(simit->second);

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
Control::init(void) throw (InitFailedException)
{
  try
  {
    _database = new Database();
    Material::set_database(*_database);

    // create the device, simulations and models
    create_device();
    create_materials();
    setup_models();
    
    // initialize the device
    _device->init();

    // initialize the simulation environments
    EnvironmentMap::iterator envit(_simulation_environments.begin());
    const EnvironmentMap::iterator envend(_simulation_environments.end());
    for ( ; envit != envend; ++envit)
      envit->second->init();


    // initialize the simulations
    SimulationMap::iterator simit(_simulations.begin());
    const SimulationMap::iterator simend(_simulations.end());
    for ( ; simit != simend; ++simit)
      simit->second->init();
  }
  catch (runtime_error& e)
  {
    string msg("Control::init failed:\nCause: ");
    msg += e.what();
    throw InitFailedException(msg);
  }
}



void
Control::create_device(void)
{
  InputParser parser(_inputfile);

  parser.read_parameters("Simulation");
  map<const string, string>& simulation_params = parser.get_parameters_map();

  if (simulation_params["searchpath"] != "")
    _database->set_search_path(simulation_params["searchpath"]);

  ModelOptions opts;
  opts["meshfile"] = simulation_params["meshfile"];
  opts["dimension"] = simulation_params["dimension"];

  _device = Device::create(opts);

  // tell the device who controls it
  _device->set_control(this);


  if (simulation_params.find("temperature") != simulation_params.end())
    SimulationOptions::temperature =
      Utils::convert<double>(simulation_params["temperature"]);
  

  cout << "Control: simulation temperature = " <<
    SimulationOptions::temperature << " K" << endl << endl;

}



void
Control::create_materials(void)
{

  assert(_device != NULL);

  InputParser parser(_inputfile);

  parser.read_device();

  map<ID, map<const string, string> >& device_map = parser.get_device_map();

  // iterate the map and create the materials
  map<ID,  map<const string, string> >::iterator mapit(device_map.begin());
  const map<ID,  map<const string, string> >::iterator mapend(device_map.end());
  for ( ; mapit != mapend; ++mapit)
  {
    map<const string, string>& data = mapit->second;
    ID region_id = mapit->first;
    const std::string& material = data["mat"];
    Material* mat = Material::create(material, ModelOptions(data));
    _device->set_material(mat, region_id);
  }
}



void
Control::setup_models(void) throw (ModelErrorException)
{

  assert(_device != NULL);

  // a reference to the device
  Device& device = *_device;
  
  InputParser parser(_inputfile);

  // parse the models section
  parser.read_models();
  map<const string, ModelStructure*>& models = parser.get_model_structure_map();

  // we loop over all simulations to setup the models
  map<const string, ModelStructure*>::iterator modit(models.begin());
  const map<const string, ModelStructure*>::iterator modend(models.end());
  for ( ; modit != modend; ++modit)
  {
    const string& modelname = modit->first;
    ModelStructure* model_str = modit->second;

    // extract the physical regions
    vector<string> preg = model_str->get_physical_regions();
    unsigned int n = preg.size();
    set<ID> phys_regions;
    for (unsigned int i = 0; i < n; i++)
      phys_regions.insert(Utils::convert<ID>(preg[i]));

    
    // read options for this simulation (from Solver and Physics section)
    parser.read_parameters("Solver", modelname);
    map<const string, string>& solveropts = parser.get_parameters_map();
    ModelOptions opts(solveropts);


    SimulationInterface* sim = SimulationInterface::create(modelname, opts);
    if (sim == NULL)
      throw ModelErrorException(
          "Control: No such simulation type: " + modelname);

    _simulations[sim->get_name()] = sim;

    // create the environment
    SimulationEnvironment* env = new SimulationEnvironment(device, phys_regions);
    _simulation_environments[sim] = env;
    sim->set_environment(env);


    //
    // now we have to create the models
    //

    opts.clear();
    // what main model should be used?
    parser.read_parameters("Physics", modelname);
    map<const string, string>& physopts = parser.get_parameters_map();
    opts += physopts;

    // parse the submodels
    map<const string, string>& physmodels = model_str->get_phys_model_map();
    opts += physmodels;

    set<ID>::iterator it(phys_regions.begin());
    const set<ID>::iterator end(phys_regions.end());

    // we have to do this for each material!
    for ( ; it != end; ++it)
    {
      Material* mat = device.get_material(*it);

      // here we actually create the model
      PhysicalModel* model = sim->create_physical_model(opts);

      // NOTE: model could be NULL, but we don't care about. Who tells us that
      // every simulation necessarily needs a model?
      mat->add_model(model, sim->get_id());
    }


    //
    // and now... the boundary conditions
    //

    map<ID, map<const string, string> >& bc_map = model_str->get_model_BC_map();
    map<ID, map<const string, string> >::iterator bdit(bc_map.begin());
    const map<ID, map<const string, string> >::iterator bdend(bc_map.end());

    for ( ; bdit != bdend; ++bdit)
    {
      ID id = bdit->first;

      opts.clear();
      opts += bdit->second;

      Boundary* bd = new Boundary(opts.get_option("BC_region_name", ""));
      BoundaryProperties* bdprop = sim->create_boundary_model(opts);

      // NOTE: bdprop could be NULL, but we don't care about. Who tells us that
      // every simulation necessarily needs a boundary model?
      bd->add_boundary_properties(bdprop, sim->get_id());
      env->add_boundary(bd, id);
    }

  } // end loop over simulations

  //
  // check for some special simulations (sweep, selfconsistency)
  // 

  // first selfconsistency

  // then sweep
  parser.read_parameters("Solver", "sweep");
  map<const string, string>& solveropts = parser.get_parameters_map();
  ModelOptions opts(solveropts);
  if (!opts.is_empty())
  {
    SimulationInterface* sim = SimulationInterface::create("sweep", opts);
    _simulations[sim->get_name()] = sim;
  }
}




void
Control::run_simulation(void) throw (SolveFailedException)
{
  InputParser parser(_inputfile);

  parser.read_parameters("Simulation");
  map<const string, string>& simulation_params = parser.get_parameters_map();

  const string& name = simulation_params["solve"];

  SimulationInterface* sim = find_simulation(name);

  if (sim == NULL)
    throw SolveFailedException("Control: Simulation not found: " + name);

  // now run it
  sim->solve();
  
}




SimulationInterface*
Control::find_simulation(const string& name) const
{
  SimulationInterface* sim = SimulationInterface::find_simulation(name);

  if (sim != NULL)
  {
    SimulationMap::const_iterator it = _simulations.find(sim->get_name());
    if (it == _simulations.end())
      sim == NULL;
  }
    
  return sim;
}



