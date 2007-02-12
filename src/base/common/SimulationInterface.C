// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Control.h"

#include "DriftDiffusion.h"
#include "ExcitonTransport.h"
#include "Macrostrain.h"
#include "Sweep.h"

#include "Utils.h"

#include "GMVIO_cell.h"

// LibMesh includes
#include "system.h"

#include <iostream>


SimulationInterface::SimulationMap
SimulationInterface::_simulation_map;



SimulationInterface::SimulationInterface(void)
  : _environment(0),
    _is_initialized(false),
    _is_solved(false),
    _relaxation_factor(1.0)
{
  ID new_id = _simulation_map.size() + 1;
  _id = new_id;

  _simulation_map[new_id] = this;
}




SimulationInterface*
SimulationInterface::create(const std::string& type,
        const ModelOptions& options)
{
  SimulationInterface* sim = NULL;

  if (type == "driftdiffusion")
    sim = DriftDiffusion::create();
  else if (type == "excitontransport")
    sim = ExcitonTransport::create();
  else if (type == "macrostrain")
    sim = Macrostrain::create();
  else if (type == "sweep")
    sim = Sweep::create();

  if (sim != NULL)
  {
    sim->set_options(options);

    // we let it know what's its identifier
    sim->set_type(type);

    //! set the name
    std::string defaultname = Utils::extract_typename(typeid(*sim));
    sim->_name = sim->get_options().get_option("name", defaultname);
    sim->_options.delete_option("name");

    std::cout << "Added simulator (ID = " << sim->get_id() <<
      " name = " << sim->get_name() << ")\n";
  }

  return sim;
}




void
SimulationInterface::destroy(SimulationInterface* p)
{
  // TODO better call a module internal method
  delete p;
}




void
SimulationInterface::init(void) throw (InitFailedException)
{
  if (!_is_initialized)
  {
    // build name for equation systems
    create_equation_system_name();
    
    do_init();
  }

  _is_initialized = true;
}




void
SimulationInterface::create_equation_system_name(void)
{
  std::ostringstream o;
  o << get_name() << get_id();
  _eq_system_name = o.str();
}




SimulationInterface*
SimulationInterface::find_simulation(const std::string& name)
{
  SimulationInterface* sim = NULL;

  SimulationMap::iterator it(_simulation_map.begin());
  SimulationMap::iterator end(_simulation_map.end());

  if (it != end)
  {
    if (name == "") // we just take the first we can find ...
      sim = it->second;
    else
    {
      // look for user defined names
      for ( ; (it != end) && ((it->second)->get_name() != name); ++it);

      if (it != end)
        sim = it->second;

      if (it == end)
      {
        // name could be model identifier
        it = _simulation_map.begin();
        for ( ; (it != end) && ((it->second)->get_type() != name); ++it);

        if (it != end)
          sim = it->second;
      }
    }
  }

  return sim;
}




EquationSystems&
SimulationInterface::get_equation_systems(void) const
{
  return _environment->get_device().get_equation_systems();
}




void
SimulationInterface::solve(void) throw (SolveFailedException) 
{
  PerfLog perflog(get_name() + ": solve", false);
  perflog.start_event("solve");

  assert(_is_initialized);

  do_solve();

  _is_solved = true;
  
  perflog.stop_event("solve");
}




NumericVector<Real>&
SimulationInterface::get_solution_vector(void)
{
  const EquationSystems& eq = get_equation_systems();
  const System& sys = eq.get_system(get_equation_system_name());

  return *sys.solution;
}



BoundaryProperties*
SimulationInterface::create_boundary_model(const ModelOptions& options) const
      throw (ModelErrorException)
{
  ignore_unused_variable(options);

  return NULL;
}


      
PhysicalModel*
SimulationInterface::create_physical_model(const ModelOptions& options) const
      throw (ModelErrorException)
{
  ignore_unused_variable(options);
  
  return NULL;
}


Control&
SimulationInterface::get_control(void)
{
  return get_environment().get_device().get_control();
}


void
SimulationInterface::plot(void)
{
  const Device& dev = get_environment().get_device();

  std::vector<double> results;
  std::vector<std::string> names;
  build_nodal_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    std::string filename = get_control().get_output_dir() + "/" +
      get_name() + "_nodal" + get_control().get_filename_suffix() + ".gmv";

    GMVIO(dev.get_mesh()).write_nodal_data(filename, results, names);
  }

  build_elemental_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    std::string filename = get_control().get_output_dir() + "/" +
      get_name() + "_elemental" + get_control().get_filename_suffix() + ".gmv";

    GMVIO_cell(dev.get_mesh()).write_ascii_cell_data(filename, results, names);
  }

}
