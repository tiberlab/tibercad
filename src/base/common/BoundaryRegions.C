// $Id$

#include "BoundaryRegions.h"

#include "elem.h"

#include <algorithm>


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
BoundaryRegions::get_bc_node_map(std::map<unsigned int, std::vector<ID> >& nodemap) const
{
  ElemSideMap::const_iterator it(_sides.begin());
  for ( ; it != _sides.end(); ++it)
  {
    ID id = it->second;
    const Elem* elem = (it->first).elem();
    unsigned int s = (it->first).side();
    AutoPtr<Elem> side(elem->build_side(s));

    for (size_t i = 0; i < side->n_nodes(); i++)
    {
      unsigned int s_id = side->node(i);
      if (std::find(nodemap[id].begin(), nodemap[id].end(), s_id) == nodemap[id].end())
        nodemap[id].push_back(s_id);
    }
  }
}
