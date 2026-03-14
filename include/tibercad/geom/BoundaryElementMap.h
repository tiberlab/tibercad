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
 * \file BoundaryElementMap.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_BOUNDARYELEMENTMAP_H
#define TC_BOUNDARYELEMENTMAP_H

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/libMeshDefs.h"

#include <map>
#include <set>
#include <vector>
#include <cstdlib> // for NULL
#include "elem.h"

class Boundary;



//! Contains a map connecting boundaries to set of elems
class BoundaryElementMap
{

  public:

    //! A typedef for convenience
    typedef std::set<const Elem*> SetType;

    //! A typedef for the map
    typedef std::map<const Boundary*, SetType> MapType;

    //! an iterator for the map
    typedef MapType::iterator map_iterator;


    //! An iterator that iterates over all elements
    /*!
     * This iterator will iterate over the elements of one or
     * all boundaries. It is possible that an element is returned
     * more than once. This happens if different contacts touch in a
     * set of nodes.
     */
    class iterator
    {
      public:

        //! Default constructor
        iterator(void) : _value(NULL) {};

        //! Constructor
        iterator(const BoundaryElementMap* container,
            const Boundary* bd = NULL, bool begin = false) :
          _container(container), _bd(bd), _value(NULL)
        {
          if (begin)
          {
            if (bd != NULL)
              _map_iter = _container->_map.find(bd);
            else
              _map_iter = _container->_map.begin();

            if (_map_iter != _container->_map.end())
            {
              _set_iter = (_map_iter->second).begin();
              if (_set_iter != (_map_iter->second).end())
                _value = &(*_set_iter);
            }
          }
        }

        //! Copy constructor
        iterator(const iterator& it) :
          _container(it._container), _bd(it._bd), _value(it._value),
          _set_iter(it._set_iter), _map_iter(it._map_iter)
        {}

        //! Prefix increment
        iterator& operator++(void)
        {
          if (_value != NULL)
          {
            _value = NULL;

            if (++_set_iter != (_map_iter->second).end())
              _value = &(*_set_iter);
            else if (_bd == NULL)
            {
              if (++_map_iter != _container->_map.end())
                if ((_set_iter = (_map_iter->second).begin()) != 
                    (_map_iter->second).end())
                  _value = &(*_set_iter);
            }
          }
          return *this;
        }

        iterator& operator=(const iterator& rhs)
        {
          _container = rhs._container;
          _bd = rhs._bd;
          _value = rhs._value;
          _set_iter = rhs._set_iter;
          _set_iter = rhs._set_iter;
          _map_iter = rhs._map_iter;
          return *this;
        }

        bool operator==(const iterator& rhs)
        {
          return (_value == rhs._value);
        }

        bool operator!=(const iterator& rhs)
        {
          return (_value != rhs._value);
        }

        const Elem* operator*(void)
        {
          return *_value;
        }

        const Boundary* get_boundary(void)
        {
          return (_value == NULL) ? NULL : _map_iter->first;
        }


      private:

        const BoundaryElementMap* _container;
        const Boundary* _bd;
        const Elem* const* _value;

        SetType::const_iterator _set_iter;
        MapType::const_iterator _map_iter;
    };





    //! Constructor
    BoundaryElementMap(void);

    //! Get the set of elements for a certain boundary 
    const SetType& get(const Boundary* boundary) const;


    //! Add a single element
    void add(const Boundary* boundary, const Elem* elem);

    //! Find the boundaries corresponding to an element
    /*!
     * \return true if the elem is on a boundary
     */
    bool find(const Elem* elem, std::set<const Boundary*>& bds) const;

    //! Find the boundaries corresponding to an element
    /*!
     * \return true if the elem is on a boundary
     */
    bool find(const Elem* elem, std::vector<const Boundary*>& bds) const;

    //! The start of the map
    map_iterator map_begin(void);

    //! The end of the map
    map_iterator map_end(void);

    //! The first element
    const iterator elements_begin(const Boundary* bd = NULL) const;

    //! The last element
    const iterator elements_end(const Boundary* bd = NULL) const;


  private:

    //! The real map
    MapType _map;

    //! The empty element set
    const static SetType _empty_set;

    //! Let the iterator access the our map
    friend class iterator;
};


//
// inline methods
//

inline
void
BoundaryElementMap::add(const Boundary* boundary, const Elem* elem)
{
  _map[boundary].insert(elem);
}


inline
BoundaryElementMap::map_iterator
BoundaryElementMap::map_begin(void)
{
  return _map.begin();
}



inline
BoundaryElementMap::map_iterator
BoundaryElementMap::map_end(void)
{
  return _map.end();
}


inline
const BoundaryElementMap::iterator
BoundaryElementMap::elements_begin(const Boundary* bd) const
{
  return iterator(this, bd, true);
}


inline
const BoundaryElementMap::iterator
BoundaryElementMap::elements_end(const Boundary* bd) const
{
  return iterator(this, bd, false);
}


#endif // TC_BOUNDARYELEMENTMAP_H
