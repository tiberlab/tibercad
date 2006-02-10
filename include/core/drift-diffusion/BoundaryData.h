// $Id$

#ifndef _BOUNDARYDATA_H_
#define _BOUNDARYDATA_H_

#include <map>
#include <vector>

// forward declarations
class Elem;
class Node;
class BoundaryDescriptor;

class BoundaryData
{
  
  public:

    typedef std::pair<const Elem*, unsigned int> ElementSide;

    typedef std::map<const ElementSide,
            BoundaryDescriptor*>::iterator iterator;

    typedef std::map<const ElementSide,
            BoundaryDescriptor*>::const_iterator const_iterator;


    int get_size(void) const;

    const_iterator sides_end(void) const;

    const_iterator sides_begin(void) const;

    const_iterator find(const ElementSide& side) const;
  
  
    void set_data(const ElementSide& side,
                  BoundaryDescriptor* descriptor);
                  
    const BoundaryDescriptor* get_data(const ElementSide& side) const;
                  
    BoundaryDescriptor* get_data(const ElementSide& side);

    BoundaryDescriptor* operator[](const ElementSide& side);

    std::vector<int> find_element(const Elem* elem) const;

    
  private:
  
    std::map<const ElementSide, BoundaryDescriptor*> _data_sides;
    
};

//
// inline member functions
//

inline
void
BoundaryData::set_data(const ElementSide& side,
                      BoundaryDescriptor* descriptor)
{
  _data_sides[side] = descriptor;
}

inline
const BoundaryDescriptor*
BoundaryData::get_data(const ElementSide& side) const
{
  const_iterator i = find(side);
  if (i != sides_end())
    return i->second;
  return 0;
}

inline
BoundaryDescriptor*
BoundaryData::get_data(const ElementSide& side)
{
  const_iterator i = find(side);
  if (i != sides_end())
    return i->second;
  return 0;
}

inline
BoundaryDescriptor*
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
