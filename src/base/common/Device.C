// $Id$

#include "Device.h"
#include "Material.h"
#include "MeshUtils.h"
#include "MeshInput.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"

#include "mesh.h"
#include "mesh_data_elements.h"
#include "equation_systems.h"

#include <iostream>

using namespace std;


Device::Device(void)
  : _mesh(NULL),
    _mesh_units(1e-6),
    _eq_system(NULL),
    _boundary_nodes(NULL),
    _symmetry(TiberCad::NONE)
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

  MeshInput::read_mesh(meshfile, dim, *_mesh, meshdata, *_boundary_nodes,
      _region_names, _boundary_region_names);

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

  const string& sym = _options.get_option("symmetry", "");
  if (sym == "cylindrical")
  {
    _symmetry = TiberCad::CYLINDRICAL;
    cout << "Using cylinder symmetry." << endl;
  }

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
    ostringstream s;
    s << "Device: physical region " << region_id <<
      " does not exist in mesh file.";
    throw InitFailedException(s.str());
  }
  if (_material_map.find(region_id) != _material_map.end())
  {
    ostringstream s;
    s << "Device: trying to redefine physical region " << region_id << ".";
    throw InitFailedException(s.str());
  }

  _material_map[region_id] = material;
  cout << "added material " << material->get_name()
    << " for region number " << region_id
    << " (" << get_region_name(region_id) << ")\n";
}


void
Device::set_material(Material* material, const std::vector<ID>& region_ids)
{
  assert(material != NULL);

  for (unsigned int i = 0; i < region_ids.size(); ++i)
    set_material(material, region_ids[i]);
}


void 
Device::get_region_ids(const std::string& name, vector<ID>& ids) const
{
  ids.resize(0);
  map<ID, string>::const_iterator it(_region_names.begin());
  const map<ID, string>::const_iterator end(_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
      ids.push_back(it->first);
}



void
Device::get_boundary_region_ids(const std::string& name, vector<ID>& ids) const
{
  ids.resize(0);
  map<ID, string>::const_iterator it(_boundary_region_names.begin());
  const map<ID, string>::const_iterator end(_boundary_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
      ids.push_back(it->first);
}



void 
Device::set_region_name(const std::string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _region_names[ids[i]] = name;
}



void
Device::set_boundary_region_name(const std::string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _boundary_region_names[ids[i]] = name;
}




double
Device::find_temperature_for_elem(const Elem* elem) const
{
  double temp = 0.0;
  
  TemperatureMap::const_iterator it(_elem_temp.find(elem));
  TemperatureMap::const_iterator end(_elem_temp.end());

  const Elem* el = elem->parent();
  while ((el != NULL) && (it != end))
  {
    el = el->parent();
    it = _elem_temp.find(el);
  }

  if (el != NULL) // we found it!
    temp = it->second;
  else
  {
    // no parent, so check for children
    // NOTE: if we find mor than one child, we build some mean value
    vector<const Elem*> tree;
    elem->family_tree(tree, false);

    unsigned int n_children = 0;
    
    unsigned int len = tree.size();
    for (unsigned int i = 0; i < len; i++)
    {
      const Elem* elem_i = tree[i];
      it = _elem_temp.find(elem_i);
      if (it != end)
      {
        temp += it->second;
        n_children++;
      }
    }

    if (n_children != 0)
      temp /= n_children;
    else
      temp = SimulationOptions::temperature;
  }

  return temp;
}



void
Device::delete_from_temperature_map(const Elem* elem)
{
  TemperatureMap::iterator end(_elem_temp.end());
  // if elem is in the list, we delete it
  while (elem != NULL)
  {
    TemperatureMap::iterator it(_elem_temp.find(elem));
    if (it != end)
    {
      _elem_temp.erase(it);
      break;
    }
    elem = elem->parent();
  }
}
