// $Id$


#include "Boundary.h"
#include "BoundaryProperties.h"
#include "SimulationEnvironment.h"

#include <cassert>


Boundary::Boundary(const std::string& name, SimulationEnvironment* environment,
    std::set<ID> region_ids)
  : _name(name),
    _area_factor(1.0),
    _env(environment)
{
  _env->add_boundary(this, region_ids);
}


void
Boundary::add_boundary_properties(BoundaryProperties* properties,
    ID simulator_id)
{
  assert(properties != NULL);
  assert(simulator_id != 0);

  // We let the model know ourself
  properties->set_boundary(this);

  properties->set_simulation_id(simulator_id);

  PropertyMap::iterator it(_models.find(simulator_id));

  if (it != _models.end())
  {
    delete it->second;
    it->second = properties;
  }
  else
    _models[simulator_id] = properties;
}

Boundary::~Boundary(void)
{
  PropertyMap::iterator it(_models.begin());
  const PropertyMap::iterator end(_models.end());

  for ( ; it != end; ++it)
    delete it->second;
}


void
Boundary::init(void)
{
  find_region_ids();

  PropertyMap::iterator it(_models.begin());
  const PropertyMap::iterator end(_models.end());

  for ( ; it != end; ++it)
    it->second->init();
  
}


void
Boundary::find_region_ids(void)
{
  assert(_env != NULL);

  SimulationEnvironment::BoundarySideIterator it(_env->boundary_sides_begin());
  const SimulationEnvironment::BoundarySideIterator end(_env->boundary_sides_end());

  for ( ; it != end; ++it)
  {
    ID id = it->second;
    if (_env->get_boundary(id) == this)
      _region_ids.insert((it->first).elem()->subdomain_id());
  }
}
