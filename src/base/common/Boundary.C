// $Id$


#include "Boundary.h"
#include "BoundaryProperties.h"

#include <cassert>


void
Boundary::add_boundary_properties(BoundaryProperties* properties,
    ID simulator_id)
{
  assert(properties != NULL);
  assert(simulator_id != 0);

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
  PropertyMap::iterator it(_models.begin());
  const PropertyMap::iterator end(_models.end());

  for ( ; it != end; ++it)
    it->second->init();
  
}
