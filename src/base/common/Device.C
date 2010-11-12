// $Id$

#include "Device.h"
#include "Material.h"
#include "MeshUtils.h"
#include "MeshReader.h"
#include "ReadGMSH.h"
#include "ReadISEGrid.h"
#include "BoundaryRegions.h"
#include "MaterialBoundary.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "MeshRegionInfo.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"

#include "mesh.h"
#include "equation_systems.h"

#include "Messages.h"

#include <iostream>

using namespace std;


namespace
{
  class IDPair
  {
    public:
      IDPair(ID a, ID b) : _a(a), _b(b) { if (b < a) {_b = a; _a = b;} }

      bool operator==(const IDPair& rhs) const { return ((_a == rhs._a) && (_b == rhs._b)); }
      bool operator!=(const IDPair& rhs) const { return !(*this == rhs); }
      bool operator<(const IDPair& rhs) const {return ((_a < rhs._a) && (_b < rhs._b)); }

    private:
      ID _a, _b;
  };
}

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

  // every ID has its own object
  BoundaryMap::iterator bdit(_boundary_map.begin());
  const BoundaryMap::iterator bdend(_boundary_map.end());
  for ( ; bdit != bdend; ++bdit)
    delete bdit->second;

  // every ID has its own object
  EdgeObjMap::iterator eit(_edge_map.begin());
  const EdgeObjMap::iterator eend(_edge_map.end());
  for ( ; eit != eend; ++eit)
    delete eit->second;

  // every ID has its own object
  NodeObjMap::iterator nit(_node_map.begin());
  const NodeObjMap::iterator nend(_node_map.end());
  for ( ; nit != nend; ++nit)
    delete nit->second;

  _material_map.clear();

  delete _eq_system;
  delete _boundary_nodes;
  delete _mesh;
  delete _mesh_region_info;
  delete _bd_regions;
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
  int dim = _options.get_option("dimension", 1);

  const string& meshfile = _options["meshfile"];

  {
    ostringstream os;
    os << "Reading mesh file " << meshfile << " ...";
    m.info(os.str(), false);
  }

  _mesh = new Mesh(dim);

  _mesh_region_info = new MeshRegionInfo;
  _bd_regions = new BoundaryRegions;

  MeshReader::read_mesh(meshfile, *_mesh, *_mesh_region_info, *_bd_regions);


  // TODO only for now (back compatibility)
  delete _boundary_nodes;
  _boundary_nodes = new map<unsigned int, vector<ID> >();
  _bd_regions->get_bc_node_map(*_boundary_nodes);


  // update mesh dimension
  dim = _mesh->mesh_dimension();

  m.info(" done.");
  {
    ostringstream os;

    os << "mesh units          : " << _mesh_units << " m" << endl
       << "mesh dimension      : " << setw(7) << setfill(' ') << dim << endl
       << "number of nodes     : " << setw(7) << setfill(' ') << _mesh->n_nodes() << endl
       << "number of elements  : " << setw(7) << setfill(' ') << _mesh->n_elem() << endl
       << "number of subdomains: " << setw(7) << setfill(' ') << _mesh_region_info->n_subdomains();
    m.info(os.str());
  }
  m.newline();

  //_mesh_region_info->print_info();

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


  _eq_system = new EquationSystems(*_mesh);
}



/*
void
Device::prepare_boundaries(void)
{
  BCNodeMap::const_iterator bd_it;
  const BCNodeMap::const_iterator bd_end(_boundary_nodes->end());

  // we only look on level zero
  MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
  const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;
    const ID id = elem->subdomain_id();

    // loop over the sides
    int n_sides = elem->n_sides();
    for (int s = 0; s < n_sides; s++)
    {

      // check if the neighbouring element is inexistent (outer boundary)
      // or in another simulation region (inner boundary)
      //
      // we allow inner 'boundaries', i.e. we don't really consider
      // boundaries but n-1 dimensional domains
      const Elem* neighbour = elem->neighbor(s);

      if ((neighbour == NULL) ||
          (neighbour->subdomain_id() != id))
      {

      }
    }
  }

}
*/


void
Device::init(void)
{

  Messages m;

  // init all materials
  m.info("Setting up bulk models ...");
  m.indent();
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it)
    (it->second)->init();
  m.unindent();
  m.newline();

  // init all lower dimensional regions
  m.info("Setting up lower dimensional (boundary) models ...");
  m.indent();
  BoundaryMap::iterator bdit(_boundary_map.begin());
  const BoundaryMap::iterator bdend(_boundary_map.end());
  for ( ; bdit != bdend; ++bdit)
    (bdit->second)->init();

  EdgeObjMap::iterator eit(_edge_map.begin());
  const EdgeObjMap::iterator eend(_edge_map.end());
  for ( ; eit != eend; ++eit)
    (eit->second)->init();

  NodeObjMap::iterator nit(_node_map.begin());
  const NodeObjMap::iterator nend(_node_map.end());
  for ( ; nit != nend; ++nit)
    (nit->second)->init();
  m.unindent();

}





void
Device::set_material(Material* material, ID region_id)
{
  assert(material != NULL);

  if (!_mesh_region_info->has_id(region_id))
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
Device::get_boundary_object(ID id)
{
  MaterialBoundary* mb = NULL;

  if (id != INVALID_ID)
  {
    BoundaryMap::iterator it(_boundary_map.find(id));
    if (it != _boundary_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_side(id))
      {
        // first get ids of the regions on both sides
        const IDSet& ids = _bd_regions->get_contiguous_regions_for_side(id);
        assert(ids.size() == 2);
        ID idA = *(ids.begin());
        ID idB = *(++(ids.begin()));
        Material* matA = get_material(idA);
        Material* matB = get_material(idB);
        mb = MaterialBoundary::create(idA, matA, idB, matB, ModelOptions());
        _boundary_map[id] = mb;
      }
    }
  }

  return mb;
}



MaterialBoundary*
Device::get_boundary_object(const Elem* elem, int side)
{
  MaterialBoundary* mb = NULL;

  ID id = _bd_regions->get_side_id(elem, side);
  if (id != INVALID_ID)
  {
    BoundaryMap::iterator it(_boundary_map.find(id));
    if (it != _boundary_map.end())
      mb = it->second;
  }

  return mb;
}


EdgeObject*
Device::get_edge_object(ID id)
{
  EdgeObject* mb = NULL;

  if (id != INVALID_ID)
  {
    EdgeObjMap::iterator it(_edge_map.find(id));
    if (it != _edge_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_edge(id))
      {
        mb = EdgeObject::create(ModelOptions());
        _edge_map[id] = mb;
      }
    }
  }

  return mb;
}


EdgeObject*
Device::get_edge_object(const Elem* elem, int edge)
{
  EdgeObject* mb = NULL;

  ID id = _bd_regions->get_edge_id(elem, edge);
  if (id != INVALID_ID)
  {
    EdgeObjMap::iterator it(_edge_map.find(id));
    if (it != _edge_map.end())
      mb = it->second;
  }

  return mb;
}



NodeObject*
Device::get_node_object(ID id)
{
  NodeObject* mb = NULL;

  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_node(id))
      {
        mb = NodeObject::create(ModelOptions());
        _node_map[id] = mb;
      }
    }
  }

  return mb;
}


NodeObject*
Device::get_node_object(const Elem* elem, int node)
{
  NodeObject* mb = NULL;

  ID id = _bd_regions->get_node_id(elem->get_node(node));
  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
  }

  return mb;
}


NodeObject*
Device::get_node_object(const Node* node)
{
  NodeObject* mb = NULL;

  ID id = _bd_regions->get_node_id(node);
  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
  }

  return mb;
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

    // next we look in the mesh region list
    if (ids.size() == 0)
    {
      const IDSet& idset = _mesh_region_info->get_ids(name);
      IDSet::iterator it(idset.begin());
      const IDSet::iterator end(idset.end());

      for ( ; it != end; ++it)
        if (_region_names.count(*it))
          ids.push_back(*it);
    }

    // as last resort we look for material names
    if (ids.size() == 0)
    {
      MaterialMap::const_iterator it(_material_map.begin());
      const MaterialMap::const_iterator end(_material_map.end());

      for ( ; it != end; ++it)
        if (it->second->get_name() == name)
          ids.push_back(it->first);
    }
  }
}



void
Device::get_mesh_region_ids(const string& name, vector<ID>& ids) const
{
  const IDSet& idset = _mesh_region_info->get_ids(name);
  ids.clear();
  ids.reserve(idset.size());

  IDSet::iterator it(idset.begin());
  const IDSet::iterator end(idset.end());
  for ( ; it != end; ++it)
    ids.push_back(*it);

  ids.reserve(ids.size());
}




void
Device::get_boundary_region_ids(const string& name, vector<ID>& ids) const
{
  IDSet idset(_bd_regions->get_ids(name));

  // if the set is empty, we try to interpret the string as a specification
  // of the boundary between two materials
  if (idset.empty())
  {
    vector<string> comp;
    Utils::tokenize(name, comp, "/");

    // there have to be two _different_ components
    if ((comp.size() == 2) && (comp[0] != comp[1]))
    {
      // it is an interface specification

      // a map to keep track of existing IDs
      map<IDPair, ID> known_ids;

      // we have to loop over all element sides!
      // we only look on level zero
      MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
      const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

      for ( ; el != el_end; ++el)
      {
        Elem* elem = *el;
        const ID id = elem->subdomain_id();

        // the material of this element
        const Material* mat = get_material(id);

        // the name of this elements region
        const string& name = get_region_name(id);

        // the other component
        size_t other;
        if ((comp[0] == mat->get_name()) || (comp[0] == name))
          other = 1;
        else if ((comp[1] == mat->get_name()) || (comp[1] == name))
          other = 0;
        else // go to the next element
          continue;


        // loop over the sides
        int n_sides = elem->n_sides();
        for (int s = 0; s < n_sides; s++)
        {

          // check if the neighbouring element is inexistent (outer boundary)
          // or in another simulation region (inner boundary)
          //
          // we allow inner 'boundaries', i.e. we don't really consider
          // boundaries but n-1 dimensional domains
          const Elem* neighbour = elem->neighbor(s);
          ID neighbour_id = INVALID_ID;

          if ((neighbour == NULL) ||
              ((neighbour_id = neighbour->subdomain_id()) != id))
          {
            map<IDPair, ID>::iterator it(known_ids.find(IDPair(id, neighbour_id)));
            // if the ID pair already exists, we can just add the elem side
            if (it != known_ids.end())
            {
              _bd_regions->add_side(elem, s, it->second);
            }
            else
            {
              if (neighbour_id == INVALID_ID)
              {

              }
              else if ((comp[other] == get_material(neighbour_id)->get_name()) ||
                       (comp[other] == get_region_name(neighbour_id)))
              {
                // first check if the side has already a boundary ID
                ID newid = _bd_regions->get_side_id(elem, s);
                if (newid == INVALID_ID)
                {
                  // we have to create a new one
                  newid = _bd_regions->next_id();
                  _bd_regions->add_side(elem, s, newid);
                  _bd_regions->set_name(newid, name);
                  known_ids[IDPair(id, neighbour_id)] = newid;
                }
                idset.insert(newid);
              }
            }
          }
        }
      }
    }
  }

  ids.clear();
  ids.reserve(idset.size());

  IDSet::iterator it(idset.begin());
  const IDSet::iterator end(idset.end());
  for ( ; it != end; ++it)
    ids.push_back(*it);
}



void
Device::set_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _region_names[ids[i]] = name;
}


/*
void
Device::set_boundary_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _boundary_region_names[ids[i]] = name;
}
*/


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



