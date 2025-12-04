/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file BoundaryNodeMap.C
 * \brief tiberCAD API implementation.
 */


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
