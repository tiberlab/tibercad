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
 * \file BoundaryElementMap.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/geom/BoundaryElementMap.h"


using namespace std;


const BoundaryElementMap::SetType
BoundaryElementMap::_empty_set;


BoundaryElementMap::BoundaryElementMap(void)
{
}


const
BoundaryElementMap::SetType&
BoundaryElementMap::get(const Boundary* boundary) const
{
  map<const Boundary*, SetType>::const_iterator it(_map.find(boundary));

  return ((it == _map.end()) ? _empty_set : it->second);
}




bool
BoundaryElementMap::find(const Elem* elem,
    set<const Boundary*>& bds) const
{
  bool found = false;

  bds.clear();

  map<const Boundary*, SetType>::const_iterator it(_map.begin());
  const map<const Boundary*, SetType>::const_iterator end(--_map.end());

  for ( ; it != end; ++it)
    if ((it->second).count(elem))
    {
      found = true;
      bds.insert(it->first);
    }

  return found;
}



bool
BoundaryElementMap::find(const Elem* elem,
    vector<const Boundary*>& bds) const
{
  bool found = false;

  bds.clear();

  map<const Boundary*, SetType>::const_iterator it(_map.begin());
  const map<const Boundary*, SetType>::const_iterator end(--_map.end());

  for ( ; it != end; ++it)
    if ((it->second).count(elem))
    {
      found = true;
      bds.push_back(it->first);
    }

  return found;
}
