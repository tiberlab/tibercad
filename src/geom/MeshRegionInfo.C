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
 * \file MeshRegionInfo.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/geom/MeshRegionInfo.h"

#include "libmesh/mesh_base.h"
#include "libmesh/parallel_implementation.h"

#include <cassert>

using namespace std;

ID
MeshRegionInfo::next_id(void) const
{
  ID maxid = 0;

  IDToNameMap::const_iterator it(_ids_to_names.begin());
  for ( ; it != _ids_to_names.end(); ++it)
    maxid = (it->first > maxid) ? it->first : maxid;

  maxid++;
  assert(maxid != INVALID_ID);

  return maxid;
}


void
MeshRegionInfo::print_info(void) const
{
  cerr << "mesh regions: " << endl;
  IDToNameMap::const_iterator it(_ids_to_names.begin());
  for ( ; it != _ids_to_names.end(); ++it)
    cerr << "  " << it->first << " : " << it->second << endl;
}


MeshRegionInfo::~MeshRegionInfo(void)
{
  clear();
}



void
MeshRegionInfo::clear(void)
{
  _ids_to_names.clear();
  _names_to_ids.clear();
}


void
MeshRegionInfo::broadcast(void)
{
  get_mesh().comm().broadcast(_ids_to_names);

  NameToIDMap::iterator it(_names_to_ids.begin());
  const NameToIDMap::iterator end(_names_to_ids.end());

  size_t mapsize = _names_to_ids.size();
  get_mesh().comm().broadcast(mapsize);

  vector<string> names(mapsize);
  if (get_mesh().comm().rank() == 0)
  {
    for (unsigned int i = 0; it != end; ++it, ++i)
      names[i] = it->first;
  }
  get_mesh().comm().broadcast(names);

  if (get_mesh().comm().rank() != 0)
  {
    for (unsigned int i = 0; i < names.size(); ++i)
      _names_to_ids[names[i]] = set<ID>();
  }

  for (it = _names_to_ids.begin(); it != end; ++it)
    get_mesh().comm().broadcast(it->second);

  do_broadcast();
}
