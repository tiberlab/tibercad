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
 * \file ElementSide.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _ELEMENTSIDE_H_
#define _ELEMENTSIDE_H_

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/libMeshDefs.h"

#include "elem.h"

#include <utility>
#include <string>



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

    ElementSide(const libMesh::Elem* elem, unsigned int side)
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


    const libMesh::Elem* elem(void) const { return _elside.first; }

    unsigned int side(void) const { return _elside.second; }


  private:

    std::pair<const libMesh::Elem*, unsigned int> _elside;
};


// we can use the same class in other contexts
typedef ElementSide ElementEdge;
typedef ElementSide ElementNode;


#endif // _ELEMENTSIDE_H_
