// $Id$

#include "Device.h"
#include "Material.h"
#include "MeshUtils.h"
#include "MeshInput.h"

#include "mesh.h"
#include "mesh_data_elements.h"
#include "equation_systems.h"

#include <iostream>

using namespace std;


Device::Device(void)
  : _mesh(NULL),
    _mesh_units(1e-6),
    _eq_system(NULL),
    _boundary_nodes(NULL)
{
  _material_map.clear();	
}


Device::~Device()
{
  // we put them first into a set because a material can be associated
  // to several IDs
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  set<Material*> mats;
  for ( ; it != end; ++it)
    mats.insert(it->second);

  set<Material*>::iterator matit(mats.begin());
  const set<Material*>::iterator matend(mats.end());
  for ( ; matit != matend; ++matit)
    delete *matit;

  _material_map.clear();

  delete _eq_system;
  delete _boundary_nodes;
  delete _mesh;
}


Device*
Device::create(const ModelOptions& options)
{

  Device* device = new Device();

  device->set_options(options);

  device->setup_mesh();

  return device;
}


void
Device::setup_mesh(void)
{

  delete _mesh;
  
  int dim = _options.get_option("dimension", 2);
  _mesh = new Mesh(dim);

  const string& meshfile = _options["meshfile"];
  
  MeshData_elements meshdata(*_mesh);
  meshdata.enable_compatibility_mode();

  delete _boundary_nodes;
  _boundary_nodes = new map<unsigned int, vector<ID> >();

  MeshInput::read_mesh(meshfile, dim, *_mesh, meshdata, *_boundary_nodes);

  MeshUtils::assign_subdomain_ids(*_mesh, meshdata);
  MeshUtils::get_subdomain_ids(*_mesh, _region_ids);

#ifdef DEBUG
  cout << endl << "Device::setup_mesh(): ";
  cout << "   meshfile : " << meshfile << " (" << dim << "D)" << endl;
  _mesh->print_info();
#endif

  _eq_system = new EquationSystems(*_mesh);
}

  
void
Device::init(void)
{
  _mesh_units = _options.get_option("mesh_units", _mesh_units);
  cout << "mesh units: " << _mesh_units << " m" << endl;

  // init all materials
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it)
    (it->second)->init();

}


void
Device::set_material(Material* material, ID region_id)
{
  assert(material != NULL);

  if (_region_ids.find(region_id) == _region_ids.end())
  {
    std::ostringstream s;
    s << "Device: physical region " << region_id <<
      " does not exist in mesh file.";
    throw InitFailedException(s.str());
  }

  _material_map[region_id] = material;
  std::cout << "added material " << material->get_name()
    << " for region number " << region_id << "\n";
}


void
Device::set_material(Material* material, const std::vector<ID>& region_ids)
{
  assert(material != NULL);

  for (unsigned int i = 0; i < region_ids.size(); ++i)
    _material_map[region_ids[i]] = material;
  
  std::cout << "added material " << material->get_name()
    << " for region numbers ";
  for (unsigned int i = 0; i < region_ids.size(); ++i)
    std::cout << " " << region_ids[i];
  std::cout << "\n";
}
