// $Id$

#include "MeshRegionInfo.h"

#include <cassert>


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
  std::cerr << "mesh regions: " << std::endl;
  IDToNameMap::const_iterator it(_ids_to_names.begin());
  for ( ; it != _ids_to_names.end(); ++it)
    std::cerr << "  " << it->first << " : " << it->second << std::endl;
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

