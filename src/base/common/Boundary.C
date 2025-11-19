// $Id$


#include "tibercad/geom/Boundary.h"
#include "tibercad/physics/BoundaryProperties.h"

#include <cassert>


Boundary::Boundary(const std::string& name, const ModelOptions& options)
  : _name(name),
    _options(options),
    _area_factor(1.0)
{
  set_area_factor(_options.get_option("area_factor", _area_factor));
}



void
Boundary::set_region_ids(const std::vector<ID>& region_ids)
{
  _region_ids.insert(region_ids.begin(), region_ids.end());
}

void
Boundary::get_region_ids(std::vector<ID>& ids) const
{
  ids.clear();
  ids.insert(ids.begin(), _region_ids.begin(), _region_ids.end());
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

  PropertyMap::iterator it(_oldmodels.find(simulator_id));

  if (it != _oldmodels.end())
  {
    delete it->second;
    it->second = properties;
  }
  else
    _oldmodels[simulator_id] = properties;
}


Boundary::~Boundary(void)
{
  PropertyMap::iterator it(_oldmodels.begin());
  const PropertyMap::iterator end(_oldmodels.end());

  for ( ; it != end; ++it)
    delete it->second;
}


void
Boundary::init(void)
{

  PropertyMap::iterator it(_oldmodels.begin());
  const PropertyMap::iterator end(_oldmodels.end());

  for ( ; it != end; ++it)
    it->second->init();

}


