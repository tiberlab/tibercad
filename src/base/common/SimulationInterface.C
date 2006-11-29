// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"

#include "DriftDiffusion.h"
#include "ExcitonTransport.h"

#include "Utils.h"

#include <iostream>


SimulationInterface::SimulationMap
SimulationInterface::_simulation_map;


SimulationInterface::SimulationInterface(void)
  : _environment(0),
    _is_initialized(false),
    _is_solved(false)
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

  if (type == "drift-diffusion")
    sim = DriftDiffusion::create();
  else if (type == "exciton-transport")
    sim = ExcitonTransport::create();


  if (sim != NULL)
  {
    sim->set_options(options);

    //! set the name
    std::string defaultname = Utils::extract_typename(typeid(*sim));
    sim->_name = sim->get_options().get_option("name", defaultname);
    sim->_options.delete_option("name");

    std::cerr << "Added simulator (ID = " << sim->get_id() <<
      " default name = " << sim->get_name() << ")\n";
  }

  return sim;
}


void
SimulationInterface::init(void) throw (InitFailedException)
{
  assert(_environment != 0);

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

  for ( ; (it != end) && ((it->second)->get_name() != name); ++it);

  if (it != end)
    sim = it->second;

  return sim;
}


EquationSystems&
SimulationInterface::get_equation_systems(void) const
{
  return _environment->get_device().get_equation_systems();
}

