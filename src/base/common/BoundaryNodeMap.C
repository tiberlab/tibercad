// $Id$

#include "tibercad/geom/BoundaryNodeMap.h"

using namespace std;



BoundaryNodeMap::BoundaryNodeMap(void)
{
  _map[INVALID_ID] = set<const libMesh::Node*>();
}


const
BoundaryNodeMap::NodeSet&
BoundaryNodeMap::get_nodes(ID id) const
{
  map<ID, NodeSet>::const_iterator it(_map.find(id));
  if (it == _map.end())
    it = _map.find(INVALID_ID);

  return it->second;
}




bool
BoundaryNodeMap::find_node(const libMesh::Node* node, set<ID>& ids) const
{
  bool found = false;

  ids.clear();

  map<ID, NodeSet>::const_iterator it(_map.begin());
  const map<ID, NodeSet>::const_iterator end(--_map.end());

  for ( ; it != end; ++it)
    if ((it->second).count(node))
    {
      found = true;
      ids.insert(it->first);
    }

  return found;
}



bool
BoundaryNodeMap::find_node(const libMesh::Node* node, vector<ID>& ids) const
{
  bool found = false;

  ids.clear();

  map<ID, NodeSet>::const_iterator it(_map.begin());
  const map<ID, NodeSet>::const_iterator end(--_map.end());

  for ( ; it != end; ++it)
    if ((it->second).count(node))
    {
      found = true;
      ids.push_back(it->first);
    }

  return found;
}
