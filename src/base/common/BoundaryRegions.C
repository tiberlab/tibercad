// $Id$

#include "BoundaryRegions.h"

#include "elem.h"

#include <algorithm>

void
BoundaryRegions::add_side(const Elem* elem, unsigned int side, ID id)
{
  _sides[ElementSide(elem, side)] = id;
  _ids_to_names[id];
}


void
BoundaryRegions::add_edge(const Elem* elem, unsigned int edge, ID id)
{
  _edges[ElementEdge(elem, edge)] = id;
  _ids_to_names[id];
}


void
BoundaryRegions::add_node(const Node* node, ID id)
{
  _nodes[node] = id;
  _ids_to_names[id];
}


void
BoundaryRegions::set_name(ID id, const std::string& name)
{
  IDToNameMap::const_iterator it(_ids_to_names.find(id));
  if (it != _ids_to_names.end())
  {
    _ids_to_names[id] = name;
    _names_to_ids[name] = id;
  }
}



ID
BoundaryRegions::get_id(const std::string& name) const
{
  ID id = INVALID_ID;
  NameToIDMap::const_iterator it(_names_to_ids.find(name));
  if (it != _names_to_ids.end())
    id = it->second;

  return id;
}



void
BoundaryRegions::clear(void)
{
  _sides.clear();
  _edges.clear();
  _nodes.clear();
  _ids_to_names.clear();
  _names_to_ids.clear();
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
