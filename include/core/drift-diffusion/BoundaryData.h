// $Id$

#ifndef _BOUNDARYDATA_H_
#define _BOUNDARYDATA_H_

#include <map>
#include <vector>

// forward declarations
class Elem;
class Node;
class ElectricalContact;

class BoundaryData
{
  
  public:

    typedef std::pair<const Elem*, unsigned int> ElementSide;

    typedef std::map<const ElementSide,
            ElectricalContact*>::iterator iterator;

    typedef std::map<const ElementSide,
            ElectricalContact*>::const_iterator const_iterator;


    int get_size(void) const;

    const_iterator sides_end(void) const;

    const_iterator sides_begin(void) const;

    const_iterator find(const ElementSide& side) const;
  
  
    void set_data(const ElementSide& side,
                  ElectricalContact* descriptor);
                  
    const ElectricalContact* get_data(const ElementSide& side) const;
                  
    ElectricalContact* get_data(const ElementSide& side);

    ElectricalContact* operator[](const ElementSide& side);

    std::vector<int> find_element(const Elem* elem) const;

    
  private:
  
    std::map<const ElementSide, ElectricalContact*> _data_sides;
    
};

//
// inline member functions
//

inline
void
BoundaryData::set_data(const ElementSide& side,
                      ElectricalContact* descriptor)
{
  _data_sides[side] = descriptor;
}

inline
const ElectricalContact*
BoundaryData::get_data(const ElementSide& side) const
{
  const_iterator i = find(side);
  if (i != sides_end())
    return i->second;
  return 0;
}

inline
ElectricalContact*
BoundaryData::get_data(const ElementSide& side)
{
  const_iterator i = find(side);
  if (i != sides_end())
    return i->second;
  return 0;
}

inline
ElectricalContact*
BoundaryData::operator[](const ElementSide& side)
{
  return _data_sides[side];
}

inline
BoundaryData::const_iterator
BoundaryData::find(const ElementSide& side) const
{
  return _data_sides.find(side);
}

inline
int
BoundaryData::get_size(void) const
{
  return _data_sides.size();
}

inline
BoundaryData::const_iterator
BoundaryData::sides_end(void) const
{
  return _data_sides.end();
}

inline
BoundaryData::const_iterator
BoundaryData::sides_begin(void) const
{
  return _data_sides.begin();
}



#endif //_BOUNDARYDATA_H_
