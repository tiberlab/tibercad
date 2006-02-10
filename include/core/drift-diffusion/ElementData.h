// $Id$

#ifndef _ELEMENTDATA_H_
#define _ELEMENTDATA_H_

#include <map>

// forward declarations
class Elem;
class SemiconductorModel;


class ElementData
{
  
  public:

    typedef std::map<const Elem*,
            SemiconductorModel*>::iterator iterator;

    typedef std::map<const Elem*,
            SemiconductorModel*>::const_iterator const_iterator;

    const_iterator end(void) const;

    const_iterator begin(void) const;

    const_iterator find(const Elem* element) const;
  
    void set_data(const Elem* element,
                  SemiconductorModel* descriptor);
                  
    SemiconductorModel* get_data(const Elem* element);

    int get_size(void) const;
    
  private:
  
    std::map<const Elem*, SemiconductorModel*> _data;
    
};

inline
void
ElementData::set_data(const Elem* element,
                      SemiconductorModel* descriptor)
{
  _data[element] = descriptor;
}

inline
SemiconductorModel*
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
