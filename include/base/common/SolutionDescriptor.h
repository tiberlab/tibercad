// $Id$


#ifndef _SOLUTIONDESCRIPTOR_H_
#define _SOLUTIONDESCRIPTOR_H_

#include "TypeDefs.h"

#include <string>
#include <cassert>
#include <iostream>


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
      TENSOR,   //!< a real 3D symmetric tensor of second rank (xx, yy, zz, xy, yz, zz)
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
     * \param units the physical units (optional)
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

    //! Get the number of components
    unsigned int n_components(void) const { return _n_comp; }

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
      assert(_n_comp > 0);
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
