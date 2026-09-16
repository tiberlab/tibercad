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
 * \file SolutionDescriptor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_SOLUTIONDESCRIPTOR_H
#define TC_SOLUTIONDESCRIPTOR_H

#include "tibercad/base/TypeDefs.h"

#include <string>
#include <cassert>
#include <iostream>
#include <set>


//! A structure describing the properties of a solution
struct SolutionDescriptor
{

  public:

    //! The numeric type of the solution variable
    enum Type
    {
      REAL,     //!< a real value
      COMPLEX,  //!< a complex value, ordered as (real, imag)
      VECTOR,   //!< a real 3D vector (x, y, z)
      TENSOR,   //!< a real 3D symmetric tensor of second rank (xx, yy, zz, xy, yz, xz)
      NTUPLE    //!< a generic n-tuple
    };


    //! The spatial association
    enum Location
    {
        NODES,  //!< located on the element nodes
        CELL,   //!< located on the cell
        ATOM,   //!< associated to an atom
        GLOBAL  //!< no spatial association
    };

    //! Default constructor
    /*!
     * Sets the ID to \c INVALID_ID
     */
    SolutionDescriptor(void) : _id(INVALID_ID) {}

    //! Constructor
    /*!
     * \param name the name used by other modules to identify a solution
     * \param id the ID to be assigned
     * \param type the numeric type
     * \param location the spatial association
     * \param units the physical units 
     * \param n_components the number of components (needed for type NTUPLE)
     */
    SolutionDescriptor(const std::string& name, ID id,
        Type type, Location location, const std::string& units = "",
        unsigned int n_components = 0);


    //! Get the ID
    /*!
     * The IDs are \em not assumed to be globally unique.
     */
    ID id(void) const { return _id; }

    //! Get the name
    const std::string& name(void) const { return _name; }

    //! Get the type
    Type type(void) const { return _type; }

    //! Get the location
    Location location(void) const { return _location; }

    //! Get the units
    const std::string& units(void) const { return _units; }

    //! Set/get the units
    std::string& units(void) { return _units; }

    //! Get the number of components
    unsigned int n_components(void) const { return _n_comp; }

    //! Set/get the number of components
    unsigned int& n_components(void) { return _n_comp; }

    //! Extract component indices from a variable name or string
    /*!
     * The syntax is [name:]xx:yy:zz, where component names or numbers are
     * divided by colons
     */
    void get_components(const std::string& str, std::set<int> comp) const;

    //! Return \c true if the quantity is located on the mesh
    bool on_mesh(void) const;

    //! To be used as map/set key
    bool operator<(const SolutionDescriptor& rhs) const
        { return this->id() < rhs.id(); };


  private:

    //! The ID assigned
    ID _id;

    //! The name
    std::string _name;

    //! The type
    Type _type;

    //! The location of the quantity inside an element
    Location _location;

    //! The units
    std::string _units;

    //! The number of components
    unsigned int _n_comp;


};


inline
SolutionDescriptor::SolutionDescriptor(const std::string& name, ID id,
    Type type, Location location, const std::string& units,
    unsigned int num_components) :
    _id(id),
    _name(name),
    _type(type),
    _location(location),
    _units(units),
    _n_comp(num_components)
{
  switch (_type)
  {
    case REAL:
      _n_comp = 1;
      break;

    case COMPLEX:
      _n_comp = 2;
      break;

    case VECTOR:
      _n_comp = 3;
      break;

    case TENSOR:
      _n_comp = 6;
      break;

    case NTUPLE:
    default:
      //assert(_n_comp > 0);
      break;
  }
}

inline
bool
SolutionDescriptor::on_mesh(void) const
{
  return ((_location == NODES) || (_location == CELL));
}


std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Type type);
std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Location location);


#endif /* _SOLUTIONDESCRIPTOR_H_ */
