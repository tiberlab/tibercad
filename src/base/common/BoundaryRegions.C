// $Id$

#include "BoundaryRegions.h"
#include "TiberCad.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "mesh_base.h"

#include <algorithm>

using namespace std;

void
BoundaryRegions::add_side(const Elem* elem, unsigned int side, ID id)
{
  if (_sides.find(ElementSide(elem, side)) == _sides.end())
  {
    _sides[ElementSide(elem, side)] = id;
    _side_ids.insert(id);
    add_id(id);

    IDSet cont_ids;
    cont_ids.insert(elem->subdomain_id());
    const Elem* neighbor = elem->neighbor(side);
    if (neighbor != NULL)
      cont_ids.insert(neighbor->subdomain_id());
    else
      cont_ids.insert(INVALID_ID);

    if (_contiguous_regions.find(id) == _contiguous_regions.end())
      _contiguous_regions[id] = cont_ids;
  }
}



void
BoundaryRegions::prepare_for_use(void)
{
  // we need a temporary map to know which ID pairs already exist
  typedef std::vector<std::set<ID> > Pairs;
  typedef HashMap<size_t, ID>::Type KnownPairs;
  typedef HashMap<ID, KnownPairs>::Type IDToKnownPairs;
  IDToKnownPairs pairs;
  Pairs pairvec;

  // We have to start from an empty _contiguous_regions structure
  _contiguous_regions.clear();

  ElemSideMap::iterator it(_sides.begin());
  const ElemSideMap::iterator end(_sides.end());
  for ( ; it != end; ++it)
  {
    ID id = it->second;
    const Elem* elem = it->first.elem();
    unsigned int side = it->first.side();

    IDSet cont_ids;
    cont_ids.insert(elem->subdomain_id());
    const Elem* neighbor = elem->neighbor(side);
    if (neighbor != NULL)
      cont_ids.insert(neighbor->subdomain_id());
    else
      cont_ids.insert(INVALID_ID);

    // find the region pair in the vector
    size_t n_pairs = pairvec.size();
    size_t pos = std::distance(pairvec.begin(),
        std::find(pairvec.begin(), pairvec.end(), cont_ids));
    if (pos == n_pairs)
      pairvec.push_back(cont_ids);

    if (_contiguous_regions.find(id) == _contiguous_regions.end())
    {
      // just add it
      _contiguous_regions[id] = cont_ids;

      pairs[id][pos] = id;
    }
    else
    {
      KnownPairs::iterator pit(pairs[id].find(pos));
      const KnownPairs::iterator pend(pairs[id].end());
      if (pit == pend)
      {
        // the region pair does not yet exist in the map
        // we have to create a new ID
        ID newid = MeshRegionInfo::next_id();
        add_id(newid);
        _side_ids.insert(newid);
        set_name(newid, get_name(id));

        it->second = newid;
        _contiguous_regions[newid] =  cont_ids;
        pairs[id][pos] = newid;
      }
      else
      {
        // the region pair already exists, just assign the right ID
        it->second = pit->second;
      }
    }
  }
}



void
BoundaryRegions::add_edge(const Elem* elem, unsigned int edge, ID id)
{
  if (_edges.find(ElementEdge(elem, edge)) == _edges.end())
  {
    _edges[ElementEdge(elem, edge)] = id;
    _edge_ids.insert(id);
    add_id(id);
  }
}


void
BoundaryRegions::add_node(const Node* node, ID id)
{
  if (_nodes.find(node) == _nodes.end())
  {
    _nodes[node] = id;
    _node_ids.insert(id);
    add_id(id);
  }
}



void
BoundaryRegions::clear(void)
{
  _sides.clear();
  _edges.clear();
  _nodes.clear();
}


void
BoundaryRegions::get_bc_node_map(std::map<ID, std::vector<ID> >& nodemap) const
{
  ElemSideMap::const_iterator it(_sides.begin());
  for ( ; it != _sides.end(); ++it)
  {
    ID id = it->second;
    const Elem* elem = (it->first).elem();
    unsigned int s = (it->first).side();
    libMesh::UniquePtr<Elem> side(elem->build_side(s));

    for (size_t i = 0; i < side->n_nodes(); i++)
    {
      unsigned int s_id = side->node(i);
      if (std::find(nodemap[id].begin(), nodemap[id].end(), s_id) == nodemap[id].end())
        nodemap[id].push_back(s_id);
    }
  }
}


/*
void
BoundaryRegions::do_broadcast(void)
{

  if (get_mesh().comm().size() == 1)
    return;

  // NOTE: TODO contiguous regions must probably be calculated globally
  // in prepare_for_use, putting together pieces from different nodes
  BDMatMap::iterator it(_contiguous_regions.begin());
  const BDMatMap::iterator end(_contiguous_regions.end());
  
  size_t reg_size = _contiguous_regions.size();
  get_mesh().comm().broadcast(reg_size);

  vector<ID> ids(reg_size);
  if (get_mesh().comm().rank() == 0)
  {
    for (unsigned int i = 0; it != end; ++it, ++i)
      ids[i] = it->first;
  }
  get_mesh().comm().broadcast(ids);
  

  if (get_mesh().comm().rank() != 0)
  {
    for (unsigned int i = 0; i < ids.size(); ++i)
      _contiguous_regions[ids[i]] = set<ID>();
  }

  for (it = _contiguous_regions.begin(); it != end; ++it)
    get_mesh().comm().broadcast(it->second);


  // now we have to distribute information on the nodes/edges/sides, but we cannot
  // distribute pointer values, so we first distribute list of dof_id_types

  // first the sides
  size_t vec_size = _sides.size();
  get_mesh().comm().broadcast(vec_size);

  vector<dof_id_type> elem_ids(vec_size);
  vector<unsigned int> elem_sides(vec_size);
  ids.resize(vec_size);

  ElemSideMap::iterator elside_it(_sides.begin());
  ElemSideMap::iterator elside_end(_sides.end());
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; elside_it != elside_end; ++elside_it, ++i)
    {
      elem_ids[i] = elside_it->first.elem()->id();
      elem_sides[i] = elside_it->first.side();
      ids[i] = elside_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(elem_sides);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _sides.clear();
    _side_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_side(mesh.elem(elem_ids[i]), elem_sides[i], ids[i]);
    }
  }

  // the the edges
  vec_size = _edges.size();
  get_mesh().comm().broadcast(vec_size);

  elem_ids.resize(vec_size);
  elem_sides.resize(vec_size);
  ids.resize(vec_size);

  elside_it = _edges.begin();
  elside_end = _edges.end();
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; elside_it != elside_end; ++elside_it, ++i)
    {
      elem_ids[i] = elside_it->first.elem()->id();
      elem_sides[i] = elside_it->first.side();
      ids[i] = elside_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(elem_sides);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _edges.clear();
    _edge_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_edge(mesh.elem(elem_ids[i]), elem_sides[i], ids[i]);
    }
  }

  // the the nodes
  vec_size = _nodes.size();
  get_mesh().comm().broadcast(vec_size);

  elem_ids.resize(vec_size);
  ids.resize(vec_size);

  NodeMap::iterator node_it(_nodes.begin());
  const NodeMap::iterator node_end(_nodes.end());
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; node_it != node_end; ++node_it, ++i)
    {
      elem_ids[i] = node_it->first->id();
      ids[i] = node_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _nodes.clear();
    _node_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_node(mesh.node_ptr(elem_ids[i]), ids[i]);
    }
  }
}
*/


void
BoundaryRegions::do_broadcast(void)
{

  if (get_mesh().comm().size() == 1)
    return;

  // NOTE: TODO contiguous regions must probably be calculated globally
  // in prepare_for_use, putting together pieces from different nodes
  BDMatMap::iterator it(_contiguous_regions.begin());
  const BDMatMap::iterator end(_contiguous_regions.end());
  
  size_t reg_size = _contiguous_regions.size();
  get_mesh().comm().broadcast(reg_size);

  vector<ID> ids(reg_size);
  if (get_mesh().comm().rank() == 0)
  {
    for (unsigned int i = 0; it != end; ++it, ++i)
      ids[i] = it->first;
  }
  get_mesh().comm().broadcast(ids);
  

  if (get_mesh().comm().rank() != 0)
  {
    for (unsigned int i = 0; i < ids.size(); ++i)
      _contiguous_regions[ids[i]] = set<ID>();
  }

  for (it = _contiguous_regions.begin(); it != end; ++it)
    get_mesh().comm().broadcast(it->second);


  // now we have to distribute information on the nodes/edges/sides, but we cannot
  // distribute pointer values, so we first distribute list of dof_id_types

  // first the sides
  size_t vec_size = _sides.size();
  get_mesh().comm().broadcast(vec_size);

  vector<libMesh::dof_id_type> elem_ids(vec_size);
  vector<unsigned int> elem_sides(vec_size);
  ids.resize(vec_size);

  ElemSideMap::iterator elside_it(_sides.begin());
  ElemSideMap::iterator elside_end(_sides.end());
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; elside_it != elside_end; ++elside_it, ++i)
    {
      elem_ids[i] = elside_it->first.elem()->id();
      elem_sides[i] = elside_it->first.side();
      ids[i] = elside_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(elem_sides);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _sides.clear();
    _side_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_side(mesh.elem(elem_ids[i]), elem_sides[i], ids[i]);
    }
  }

  // the the edges
  vec_size = _edges.size();
  get_mesh().comm().broadcast(vec_size);

  elem_ids.resize(vec_size);
  elem_sides.resize(vec_size);
  ids.resize(vec_size);

  elside_it = _edges.begin();
  elside_end = _edges.end();
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; elside_it != elside_end; ++elside_it, ++i)
    {
      elem_ids[i] = elside_it->first.elem()->id();
      elem_sides[i] = elside_it->first.side();
      ids[i] = elside_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(elem_sides);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _edges.clear();
    _edge_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_edge(mesh.elem(elem_ids[i]), elem_sides[i], ids[i]);
    }
  }

  // the the nodes
  vec_size = _nodes.size();
  get_mesh().comm().broadcast(vec_size);

  elem_ids.resize(vec_size);
  ids.resize(vec_size);

  NodeMap::iterator node_it(_nodes.begin());
  const NodeMap::iterator node_end(_nodes.end());
  if (get_mesh().comm().rank() == 0)
  {
    for (size_t i = 0; node_it != node_end; ++node_it, ++i)
    {
      elem_ids[i] = node_it->first->id();
      ids[i] = node_it->second;
    }
  }
  get_mesh().comm().broadcast(elem_ids);
  get_mesh().comm().broadcast(ids);

  // now we can reconstruct the actual information
  if (get_mesh().comm().rank() != 0)
  {
    const MeshBase& mesh = get_mesh();
    _nodes.clear();
    _node_ids.clear();

    for (size_t i = 0; i < elem_ids.size(); ++i)
    {
      add_node(mesh.node_ptr(elem_ids[i]), ids[i]);
    }
  }
}
