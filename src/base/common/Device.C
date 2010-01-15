// $Id$

#include "Device.h"
#include "Material.h"
#include "MeshUtils.h"
#include "MeshReader.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"

#include "mesh.h"
#include "mesh_data_elements.h"
#include "equation_systems.h"

#include "Messages.h"

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
  delete _meshdata;
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

  Messages m;
  m.info("Setup mesh");
  m.indent();

  _mesh_units = _options.get_option("mesh_units", _mesh_units);

  // this is backup solution if dim cannot be guessed from the mesh file
  int dim = _options.get_option("dimension", -1);

  const string& meshfile = _options["meshfile"];

  {
    ostringstream os;
    os << "Reading mesh file " << meshfile << " ...";
    m.info(os.str(), false);
  }


  delete _boundary_nodes;
  _boundary_nodes = new map<unsigned int, vector<ID> >();

  MeshReader::read_mesh(meshfile, dim, _mesh, _meshdata, *_boundary_nodes,
                       _mesh_region_names, _boundary_region_names);


  MeshUtils::assign_subdomain_ids(*_mesh, *_meshdata);

  MeshUtils::get_subdomain_ids(*_mesh, _region_ids);

  // update mesh dimension
  dim = _mesh->mesh_dimension();

  m.info(" done.");
  {
    ostringstream os;

    os << "mesh units          : " << _mesh_units << " m" << endl
       << "mesh dimension      : " << setw(7) << setfill(' ') << dim << endl
       << "number of nodes     : " << setw(7) << setfill(' ') << _mesh->n_nodes() << endl
       << "number of elements  : " << setw(7) << setfill(' ') << _mesh->n_elem() << endl;
    //<< "number of subdomains: " << setw(7) << setfill(' ') << _mesh->n_subdomains();
    m.info(os.str());
  }
  m.newline();

  const string& sym = _options.get_option("symmetry", "");
  if (sym == "cylindrical")
  {
    _symmetry = TiberCad::CYLINDRICAL;
    m.info("Using cylinder symmetry (=> 3D simulation)");
  }
  

  /*
   * NOTE:
   * In parallel case, the local mesh does not contain all elements
   * and therefore only a part of the region IDs are present
   */

#ifdef DEBUG
  cout << endl << "Device::setup_mesh(): ";
  _mesh->print_info();
#endif

  _eq_system = new EquationSystems(*_mesh);
}

  



void
Device::init(void)
{

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
    /*
     * In single processor case this has to be considered an error,
     * in parallel we should do the check in a different manner
     */
    if (libMesh::n_processors() == 1)
    {
      ostringstream s;
      s << "Device: region " << region_id <<
        " does not exist in mesh file.";
      throw InitFailedException(s.str());
    }
  }
  if (_material_map.find(region_id) != _material_map.end())
  {
    ostringstream s;
    s << "Device: trying to redefine mesh region " << region_id << ".";
    throw InitFailedException(s.str());
  }

  _material_map[region_id] = material;
  _active_region_ids.insert(region_id);
}



void
Device::set_material(Material* material, const vector<ID>& region_ids,
                     const string& region_name)
{
  assert(material != NULL);

  for (unsigned int i = 0; i < region_ids.size(); ++i)
    set_material(material, region_ids[i]);

  set_region_name(region_name, region_ids);

  ostringstream os;
  os << "Added material " << material->get_name()
     << " for region \'" << region_name
     << "\' (mesh regions " << region_ids[0];
  for (unsigned int i = 1; i < region_ids.size(); ++i)
    os << ", " << region_ids[i];
  os << ")" << endl;
  Messages::info(os.str());
}





Material*
Device::get_material(const std::string& name)
{
  Material* mat = NULL;

  map<ID, string>::const_iterator it(_region_names.begin());
  const map<ID, string>::const_iterator end(_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
    {
      mat = get_material(it->first);
      break;
    }

  return mat;
}


const Material*
Device::get_material(const std::string& name) const
{
  const Material* mat = NULL;

  map<ID, string>::const_iterator it(_region_names.begin());
  const map<ID, string>::const_iterator end(_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
    {
      mat = get_material(it->first);
      break;
    }

  return mat;
}


MaterialBoundary*
Device::get_boundary_object(const Elem* elem, int side)
{
  return NULL;
}


EdgeObject*
Device::get_edge_object(const Elem* elem, int edge)
{
  return NULL;
}


NodeObject*
Device::get_node_object(const Elem* elem, int node)
{
  return NULL;
}



void 
Device::get_active_region_ids(const string& name, vector<ID>& ids) const
{
  ids.resize(0);

  ClusterMap::const_iterator clit(_cluster_map.find(name));
  if (clit != _cluster_map.end())
    ids = clit->second;
  else
  {
    map<ID, string>::const_iterator it(_region_names.begin());
    const map<ID, string>::const_iterator end(_region_names.end());
    for ( ; it != end; ++it)
      if (it->second == name)
        ids.push_back(it->first);

    // as last resort we look in the mesh region list
    if (ids.size() == 0)
    {
      map<ID, string>::const_iterator it(_mesh_region_names.begin());
      const map<ID, string>::const_iterator end(_mesh_region_names.end());
      for ( ; it != end; ++it)
        if ((it->second == name) && (_active_region_ids.count(it->first) == 1))
          ids.push_back(it->first);
    }
  }
}



void 
Device::get_mesh_region_ids(const string& name, vector<ID>& ids) const
{
  ids.resize(0);

  map<ID, string>::const_iterator it(_mesh_region_names.begin());
  const map<ID, string>::const_iterator end(_mesh_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
      ids.push_back(it->first);
}




void
Device::get_boundary_region_ids(const string& name, vector<ID>& ids) const
{
  ids.resize(0);
  map<ID, string>::const_iterator it(_boundary_region_names.begin());
  const map<ID, string>::const_iterator end(_boundary_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
      ids.push_back(it->first);
}



void 
Device::set_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _region_names[ids[i]] = name;
}



void
Device::set_boundary_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _boundary_region_names[ids[i]] = name;
}



void
Device::set_cluster(const string& name, const vector<ID>& ids)
{
  if (_cluster_map.find(name) != _cluster_map.end()) {
    string msg("Cluster ");
    msg += name;
    msg += " already defined.";
    throw InitFailedException(msg);
  }

  _cluster_map[name] = ids;
}



