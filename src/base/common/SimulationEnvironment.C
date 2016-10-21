// $Id$

#include "SimulationEnvironment.h"
#include "Control.h"
#include "Device.h"
#include "Boundary.h"

#include "mesh.h"
#include "elem.h"
#include "libmesh_logging.h"

#include <vector>


using namespace std;

//USELIBMESHTYPE(UniquePtr);


SimulationEnvironment::EnvironmentSet
SimulationEnvironment::_environments;


SimulationEnvironment::SimulationEnvironment(
    Device& device, const set<ID>& region_numbers)
  : _device(&device),
    _mesh(&device.get_mesh()),
    _region_numbers(region_numbers),
    _is_initialized(false),
    _is_prepared(false)
{
  _environments.insert(this);
}



SimulationEnvironment::~SimulationEnvironment(void)
{
  BoundaryIterator it(boundaries_begin());
  const BoundaryIterator end(boundaries_end());
  for ( ; it != end; ++it)
    delete *it;

  _environments.erase(this);
}

void
SimulationEnvironment::set_mesh(MeshBase* mesh)
{
  if (_mesh != mesh)
  {
    _mesh = mesh;
    if (this->is_prepared())
      prepare();
  }
}


void
SimulationEnvironment::prepare(void)
{
  assert(_device != NULL);
  create_element_list();
  create_bc_maps();
  update_boundary_node_map();
  calculate_bounding_box();
}



void
SimulationEnvironment::init(void)
{
  assert(_device != NULL);

  if (!_is_initialized)
  {
    BCMap::iterator it = _bc_map.begin();
    const BCMap::iterator end = _bc_map.end();

    for ( ; it != end; ++it)
      (*it)->init();

    _is_initialized = true;
  }
}



void
SimulationEnvironment::add_boundary(Boundary* boundary)
{
  if (boundary == NULL) return;

  _bc_map.insert(boundary);

  // TODO somewhere we should check if a boundary has already been assigned
}


void
SimulationEnvironment::create_element_list(void)
{
  _element_list.clear();
  MeshBase& mesh = get_mesh();

  MeshBase::element_iterator it = mesh.active_local_elements_begin();
  const MeshBase::element_iterator end = mesh.active_local_elements_end();

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
  START_LOG("create_bc_maps", "");

  const Device::BCNodeMap& bd_nodes = _device->get_boundary_node_map();
  Device::BCNodeMap::const_iterator bd_it;
  const Device::BCNodeMap::const_iterator bd_end(bd_nodes.end());

  BCMap::const_iterator bc_it;
  const BCMap::const_iterator bc_end(_bc_map.end());

  const MeshBase& mesh = get_mesh();
  const unsigned dim = mesh.mesh_dimension();

  // we only look on level zero
  MeshBase::const_element_iterator el(mesh.local_level_elements_begin(0));
  const MeshBase::const_element_iterator el_end(mesh.local_level_elements_end(0));

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
            const set<ID>& bd_ids = (*bc_it)->get_region_ids();

            set<ID>::const_iterator it(bd_ids.begin());
            for ( ; it != bd_ids.end(); ++it)
            {
              ID bd_id = *it;

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
                  libMesh::UniquePtr<Elem> side = elem->build_side(s);
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
  }

  STOP_LOG("create_bc_maps", "");
}





const SimulationEnvironment::BoundarySideIterator
SimulationEnvironment::boundary_sides_begin(const std::string& name) const
{
  if (name.empty())
    return BoundarySideIterator(_element_side_map, _element_side_map.begin());

  return BoundarySideIterator(_element_side_map, _element_side_map.begin(),
      get_boundary(name)->get_region_ids());
}




const SimulationEnvironment::BoundarySideIterator
SimulationEnvironment::boundary_sides_end(const std::string& name) const
{
  if (name.empty())
    return BoundarySideIterator(_element_side_map, _element_side_map.end());

  return BoundarySideIterator(_element_side_map, _element_side_map.end(),
      get_boundary(name)->get_region_ids());
}




void
SimulationEnvironment::update_boundary_node_map(void)
{
  // If nothing strange happens, we should never loose boundary nodes but
  // just add new ones. So we don't reset the map

  BoundarySideIterator it(boundary_sides_begin());
  const BoundarySideIterator end(boundary_sides_end());

  if ((get_mesh()).mesh_dimension() == 1)
  {
    // 1D case is easy: boundary nodes will always be the same
    // Nevertheless we always compute them as this does not take much time
    for ( ; it != end; ++it)
    {
      const ElementSide& elem_side = it->first;
      const Elem* elem = elem_side.elem();
      //_node_map[(elem_side.first)->get_node(elem_side.second)] = it->second;
      _node_map.add_node(it->second,
          (elem_side.elem())->get_node(elem_side.side()));
    }
  }
  else
  {
    for ( ; it != end; ++it)
    {
      const ElementSide& elem_side = it->first;
      const Elem* elem = elem_side.elem();
      const unsigned int side_num = elem_side.side();

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

        ElementSide elemside(child, side_num);
        // we allow for internal boundaries
        if (is_boundary(elemside))
        {
          libMesh::UniquePtr<Elem> side = child->build_side(side_num);
          for (unsigned int i = 0; i < side->n_nodes(); i++)
            //_node_map[side->get_node(i)] = it->second;
            _node_map.add_node(it->second, side->get_node(i));
        }
      }
    }
  }
}




void
SimulationEnvironment::update_boundary_element_map(
    set<const Boundary*> boundaries)
{
  if (boundaries.size() == 0)
    boundaries.insert(_bc_map.begin(), _bc_map.end());

  set<const Boundary*>::iterator ctit;
  const set<const Boundary*>::iterator ctend(boundaries.end());

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it(mesh.active_elements_begin());
  const MeshBase::element_iterator end(mesh.active_elements_end());

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    vector<ID> bd_set;
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
      for (ctit = boundaries.begin(); ctit != ctend; ++ctit)
        if (is_node_on_boundary(elem->get_node(n), *ctit))
          _bd_elem_map.add(*ctit, elem);
  }

}





Boundary*
SimulationEnvironment::get_boundary(const string& name) const
{
  Boundary* bd = NULL;

  BCMap::const_iterator it(_bc_map.begin());
  const BCMap::const_iterator end(_bc_map.end());

  for ( ; it != end; ++it)
    if ((*it)->get_name() == name)
      bd = *it;

    return bd;
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

  MeshBase& mesh = get_mesh();

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

  invalidate_all();

  _is_prepared = true;
}



void
SimulationEnvironment::invalidate_all(void)
{
  EnvironmentSet::iterator it(_environments.begin());
  const EnvironmentSet::iterator end(_environments.end());

  for ( ; it != end; ++it)
    (*it)->invalidate();
}


void
SimulationEnvironment::destroy(SimulationEnvironment* env)
{
  if (_environments.find(env) != _environments.end())
    delete env;
}


bool
SimulationEnvironment::is_node_on_boundary(const Node* node,
    const Boundary* boundary) const
{
  bool found = false;

  BCMap::const_iterator it(_bc_map.find(const_cast<Boundary*>(boundary)));
  const BCMap::const_iterator end(_bc_map.end());

  if (it != end)
  {
    const set<ID>& ids = (*it)->get_region_ids();
    set<ID>::const_iterator idit(ids.begin());
    for ( ; idit != ids.end(); ++idit)
    {
      const BoundaryNodeMap::NodeSet& nodes = _node_map.get_nodes(*idit);
      if (nodes.count(node))
      {
        found = true;
        break;
      }
    }
  }

  return found;
}



void
SimulationEnvironment::calculate_bounding_box(void)
{
  MeshBase::const_element_iterator it = get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end = get_mesh().active_local_elements_end();

  Point pmax(std::numeric_limits<double>::min(),
             std::numeric_limits<double>::min(),
             std::numeric_limits<double>::min());
  Point pmin(std::numeric_limits<double>::max(),
             std::numeric_limits<double>::max(),
             std::numeric_limits<double>::max());
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    if (!_region_numbers.count(elem->subdomain_id()))
      continue;

    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
    {
      const Point& p = elem->point(i);
      if (p(0) < pmin(0))
        pmin(0) = p(0);
      else if (p(0) > pmax(0))
        pmax(0) = p(0);

      if (p(1) < pmin(1))
        pmin(1) = p(1);
      else if (p(1) > pmax(1))
        pmax(1) = p(1);

      if (p(2) < pmin(2))
        pmin(2) = p(2);
      else if (p(2) > pmax(2))
        pmax(2) = p(2);

    }
  }

  if (get_mesh().comm().size() > 1)
  {
    vector<double> point = {pmin(0), pmin(1), pmin(2)};
    get_mesh().comm().min(point);
    pmin = Point(point[0], point[1], point[2]);

    point = {pmax(0), pmax(1), pmax(2)};
    get_mesh().comm().max(point);
    pmax = Point(point[0], point[1], point[2]);
  }

  _bbox = make_pair(pmin, pmax);
}


std::pair<Point, Point>
SimulationEnvironment::get_bounding_box(bool recalculate)
{
  if (recalculate)
    calculate_bounding_box();

  return(_bbox);
}


Boundary*
SimulationEnvironment::get_boundary(ID boundary_number) const
{
  BCMap::const_iterator it(_bc_map.begin());
  for ( ; it != _bc_map.end(); ++it)
    if ((*it)->has_region_id(boundary_number))
      return *it;

  return NULL;
}



void
SimulationEnvironment::get_boundary_nodes(const Boundary* boundary,
    set<const Node*>& nodelist)
{
  nodelist.clear();

  if (boundary)
  {
    vector<ID> ids;
    boundary->get_region_ids(ids);

    for (size_t i = 0; i < ids.size(); ++i)
    {
      const set<const Node*>& list = _node_map.get_nodes(ids[i]);
      nodelist.insert(list.begin(), list.end());
    }
  }
}


void
SimulationEnvironment::get_region_names(set<string>& names) const
{
  names.clear();
  const set<ID>& ids = get_region_ids();
  for (set<ID>::const_iterator it(ids.begin()); it != ids.end(); ++it)
  {
    names.insert(get_device().get_region_name(*it));
  }
}
