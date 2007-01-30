// $Id$

#include "Device.h"
#include "Material.h"

#include "equation_systems.h"

#include <iostream>


Device::Device(Mesh& mesh, BoundaryNodeMap& boundary_nodes,
    const ModelOptions& options)
  : _mesh(&mesh),
    _boundary_nodes(&boundary_nodes),
    _mesh_units(1e-6),
    _options(options)
{
  _material_map.clear();	
  _eq_system = new EquationSystems(mesh);
}


Device::~Device()
{
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it)
    delete it->second;

  _material_map.clear();

  delete _eq_system;
}

void
Device::init(void)
{
  _mesh_units = _options.get_option("mesh_units", _mesh_units);

  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it)
    (it->second)->init();
}


void
Device::set_material(Material* material, ID region_id)
{
  assert(material != NULL);

  _material_map[region_id] = material;
  std::cerr << "added material " << material->get_name()
    << " for region number " << region_id << "\n";
}


void
Device::set_material(Material* material, const std::vector<ID>& region_ids)
{
  assert(material != NULL);

  for (unsigned int i = 0; i < region_ids.size(); ++i)
    _material_map[region_ids[i]] = material;
  
  std::cerr << "added material " << material->get_name()
    << " for region numbers ";
  for (unsigned int i = 0; i < region_ids.size(); ++i)
    std::cerr << " " << region_ids[i];
  std::cerr << "\n";
}
