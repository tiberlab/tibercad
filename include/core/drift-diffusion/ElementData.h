// $Id$

#ifndef _ELEMENTDATA_H_
#define _ELEMENTDATA_H_

#include <map>
#include <cassert>

// forward declarations
class Elem;
class DriftDiffusionProperties;


class ElementData
{
  
  public:

    typedef std::map<const Elem*,
            DriftDiffusionProperties*>::iterator iterator;

    typedef std::map<const Elem*,
            DriftDiffusionProperties*>::const_iterator const_iterator;

    const_iterator end(void) const;

    const_iterator begin(void) const;

    const_iterator find(const Elem* element) const;
  
    void set_data(const Elem* element,
                  DriftDiffusionProperties* descriptor);
                  
    DriftDiffusionProperties* get_data(const Elem* element);

    int get_size(void) const;
    
  private:
  
    std::map<const Elem*, DriftDiffusionProperties*> _data;
    
};

inline
void
ElementData::set_data(const Elem* element,
                      DriftDiffusionProperties* descriptor)
{
  _data[element] = descriptor;
}

inline
DriftDiffusionProperties*
ElementData::get_data(const Elem* element)
{
  assert(_data.find(element) != _data.end());

  return _data[element];
}

inline
ElementData::const_iterator
ElementData::find(const Elem* element) const
{
  return _data.find(element);
}

inline
int
ElementData::get_size(void) const
{
  return _data.size();
}

inline
ElementData::const_iterator
ElementData::end(void) const
{
  return _data.end();
}

inline
ElementData::const_iterator
ElementData::begin(void) const
{
  return _data.begin();
}


#endif //_ELEMENTDATA_H_
