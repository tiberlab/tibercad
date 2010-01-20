// $Id$

#ifndef _ELEMENTSIDE_H_
#define _ELEMENTSIDE_H_

#include <utility>


class Elem;


//! A class defining an element side
class ElementSide
{

  public:

    struct hash 
    {
      size_t operator()(const ElementSide& elside) const
      {
        size_t x = reinterpret_cast<size_t>(elside._elside.first);
        unsigned int y = elside._elside.second;
        return (x << 4) | (y && 0x000f);
      }
    };

    ElementSide(const Elem* elem, unsigned int side)
      : _elside(elem, side) {}

    bool operator==(const ElementSide& rhs) const
    {
      //return (this->elem() == rhs.elem()) && (this->side() == rhs.side());
      return (this->_elside == rhs._elside);
    }

    bool operator<(const ElementSide& rhs) const
    {
      return (this->_elside < rhs._elside);
    }


    const Elem* elem(void) const { return _elside.first; }

    unsigned int side(void) const { return _elside.second; }


  private:

    std::pair<const Elem*, unsigned int> _elside;
};

#endif // _ELEMENTSIDE_H_
