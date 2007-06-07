// $Id$

#include "SimulationEnvironment.h"
#include "Control.h"
#include "Device.h"
#include "Boundary.h"

#include "mesh.h"
#include "elem.h"

#include <vector>


using namespace std;

SimulationEnvironment::SimulationEnvironment(
    Device& device, set<ID> region_numbers)
  : _device(&device),
    _region_numbers(region_numbers),
    _is_initialized(false),
    _is_prepared(false)
{
}


SimulationEnvironment::SimulationEnvironment(
    Device& device, ID region_number)
  : _device(&device),
    _is_initialized(false),
    _is_prepared(false)
{
  _region_numbers.insert(region_number);
}


SimulationEnvironment::~SimulationEnvironment(void)
{
  BCMap::iterator it = _bc_map.begin();
  const BCMap::iterator end = _bc_map.end();

  // we have to put them first into a set, because each boundary
  // can be associated to several IDs
  set<Boundary*> bds;
  for ( ; it != end; ++it)
    bds.insert(it->second);

  set<Boundary*>::iterator bdit(bds.begin());
  const set<Boundary*>::iterator bdend(bds.end());
  for ( ; bdit != bdend; ++bdit)
    delete *bdit;
}

void
SimulationEnvironment::init(void)
{
  assert(_device != NULL);

  if (!_is_initialized)
  {
    create_element_list();
    create_bc_maps();
    update_boundary_node_map();

    BCMap::iterator it = _bc_map.begin();
    const BCMap::iterator end = _bc_map.end();

    for ( ; it != end; ++it)
      it->second->init();

    _is_initialized = true;
  }
}

void
SimulationEnvironment::add_boundary(Boundary* boundary, ID boundary_id)
{
  assert(boundary != NULL);
  assert(!_is_initialized);

  if ((_device->get_boundary_node_map()).find(boundary_id) == 
      (_device->get_boundary_node_map()).end())
  {
    std::ostringstream s;
    s << "SimulationEnvironment: boundary region " << boundary_id <<
      " does not exist in device.";
    throw InitFailedException(s.str());
  }
  
  BCMap::iterator it = _bc_map.find(boundary_id);
  if (it == _bc_map.end())
    _bc_map[boundary_id] = boundary;
  else
  {
    std::ostringstream s;
    s << "SimulationEnvironment: boundary region " << boundary_id <<
      " already defined.";
    throw InitFailedException(s.str());
    //delete it->second;
    //it->second = boundary;
  }
}

void
SimulationEnvironment::add_boundary(Boundary* boundary,
    const set<ID>& boundary_ids)
{
  assert(boundary != NULL);
  
  set<ID>::const_iterator it(boundary_ids.begin());
  set<ID>::const_iterator end(boundary_ids.end());

  for ( ; it != end; ++it)
    add_boundary(boundary, *it);
}


void
SimulationEnvironment::create_element_list(void)
{
  _element_list.clear();

  Mesh& mesh = _device->get_mesh();

  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  const set<ID>::iterator list_end = _region_numbers.end();

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;

    // to be sure we check for the subdomain id
    if (_region_numbers.find(elem->subdomain_id()) != list_end)
      _element_list.insert(elem);
  }
}



void
SimulationEnvironment::create_bc_maps(void)
{
  PerfLog perf_log("create_bc_maps", false);
  perf_log.start_event("create_bc_maps");

  const BoundaryNodeMap& bd_nodes = _device->get_boundary_node_map();
  BoundaryNodeMap::const_iterator bd_it;
  const BoundaryNodeMap::const_iterator bd_end(bd_nodes.end());

  BCMap::const_iterator bc_it;
  const BCMap::const_iterator bc_end(_bc_map.end());

  const Mesh& mesh = _device->get_mesh();
  const unsigned dim = mesh.mesh_dimension();

  // we only look on level zero
  Mesh::const_element_iterator el(mesh.level_elements_begin(0));
  const Mesh::const_element_iterator el_end(mesh.level_elements_end(0));

  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;
    const ID id = elem->subdomain_id();

    // only if the element is inside of our simulation region it can have
    // a boundary side
    if (this->contains_region(id))
    {
      
      int n_sides = elem->n_sides();
      for (int s = 0; s < n_sides; s++)
      {

        // check if the neighbouring element is inexistent (outer boundary)
        // or in another simulation region (inner boundary)
        //
        // we allow inner 'boundaries', i.e. we don't really consider
        // boundaries but n-1 dimensional domains
        if (is_boundary(ElementSide(elem, s)))
        {
          
          // now we also have to loop over all relevant map entries
          // in bd_nodes. We look only for the IDs that are found in the
          // PropertyMap
          for (bc_it = _bc_map.begin(); bc_it != bc_end; ++bc_it)
          {
            const ID bd_id = bc_it->first;

            bd_it = bd_nodes.find(bd_id);
            if (bd_it != bd_end)
            {
              const vector<ID>& nodes = bd_it->second;
              const vector<ID>::const_iterator n_begin(nodes.begin());
              const vector<ID>::const_iterator n_end(nodes.end());
              
              if (dim == 1)
              {
                // we cannot build the element sides here
                if (find(n_begin, n_end, elem->node(s)) != n_end)
                  _element_side_map[ElementSide(elem, s)] = bd_id;
              }
              else
              {
                bool found = true;
                AutoPtr<Elem> side = elem->build_side(s);
                // check if all nodes of the side are in the node map
                for (unsigned int i = 0; i < side->n_nodes(); i++)
                {
                  if (find(n_begin, n_end, side->node(i)) == n_end)
                    found = false;
                }
                if (found)
                  _element_side_map[ElementSide(elem, s)] = bd_id;
              }
            }
          }
        }
      }
    }
  }
  perf_log.stop_event("create_bc_maps");
}


void
SimulationEnvironment::update_boundary_node_map(void)
{
  // If nothing strange happens, we should never loose boundary nodes but
  // just add new ones. So we don't reset the map

  BoundarySideIterator it(_element_side_map.begin());
  const BoundarySideIterator end(_element_side_map.end());
  
  if ((_device->get_mesh()).mesh_dimension() == 1)
  {
    // 1D case is easy: boundary nodes will always be the same
    // Nevertheless we always compute them as this does not take much time
    for ( ; it != end; ++it)
    {
      const ElementSide& elem_side = it->first;
      _node_map[(elem_side.first)->get_node(elem_side.second)] = it->second;
    }
  }
  else
  {
    for ( ; it != end; ++it)
    {
      const ElementSide& elem_side = it->first;
      const Elem* elem = elem_side.first;
      const unsigned int side_num = elem_side.second;

      // get the active family tree of this element
      vector<const Elem*> fam_tree;
      elem->active_family_tree(fam_tree);

      // loop over all active children and find boundary sides that
      // correspond to side_num
      vector<const Elem*>::const_iterator elem_it(fam_tree.begin());
      vector<const Elem*>::const_iterator fam_end(fam_tree.end());
      for ( ; elem_it != fam_end; ++elem_it)
      {
        const Elem* child = *elem_it;

        // we allow for internal boundaries
        if (is_boundary(ElementSide(child, side_num)))
        {
          AutoPtr<Elem> side = child->build_side(side_num);
          for (unsigned int i = 0; i < side->n_nodes(); i++)
            _node_map[side->get_node(i)] = it->second;
        }
      }
    }
  }
}


Boundary*
SimulationEnvironment::get_boundary(const string& name) const
{
  Boundary* bd = NULL;

  BCMap::const_iterator it(_bc_map.begin());
  const BCMap::const_iterator end(_bc_map.end());

  for ( ; it != end; ++it)
    if ((it->second)->get_name() == name)
      bd = it->second;

    return bd;
}


ID
SimulationEnvironment::get_boundary_id(const string& name) const
{
  ID id = 0;

  BCMap::const_iterator it(_bc_map.begin());
  const BCMap::const_iterator end(_bc_map.end());

  for ( ; it != end; ++it)
    if ((it->second)->get_name() == name)
      id = it->first;

  return id;
}


void
SimulationEnvironment::prepare_for_solve(void)
{
  assert(_device != NULL);

  // check if we are already prepared
  if (is_prepared())
    return;

  if (!_is_initialized)
    init();
  
  Mesh& mesh = _device->get_mesh();

  MeshBase::element_iterator it = mesh.elements_begin();
  const MeshBase::element_iterator end = mesh.elements_end();

  for ( ; it != end; ++it)
  {
    Elem* el = *it;

    if (contains_element(el))
      el->set_refinement_flag(Elem::DO_NOTHING);
    else
      el->set_refinement_flag(Elem::INACTIVE);
  }

  (_device->get_control()).invalidate_environments();
  
  _is_prepared = true;
}


void
SimulationEnvironment::get_boundary_nodes(const Boundary* boundary,
    set<const Node*>& nodelist)
{
  nodelist.clear();

  BoundaryNodeIterator it(boundary_nodes_begin());
  const BoundaryNodeIterator end(boundary_nodes_end());

  for ( ; it != end; ++it)
  {
    // we can assume that any ID we get here is also in the BCMap
    if (_bc_map[it->second] == boundary)
      nodelist.insert(it->first);
  }
}
