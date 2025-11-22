// $Id$

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
