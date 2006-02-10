// $Id: BoundaryData.C 4 2005-10-28 12:05:05Z maufder $

#include "BoundaryData.h"

#include "elem.h"

std::vector<int>
BoundaryData::find_element(const Elem* elem) const
{
  std::vector<int> sides;

  const_iterator it = _data_sides.begin();
  for ( ; it != _data_sides.end(); ++it)
  {
    if ((it->first).first == elem)
      sides.push_back((it->first).second);
  }


  return sides;
}
