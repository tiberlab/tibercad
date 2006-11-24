// $Id$

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "TypeDefs.h"

#include <map>
#include <string>

class BoundaryProperties;


//! The class that contains the models for a simulation boundary
/*!
 * A simulation boundary can be an external device boundary, but also
 * an internal boundary between simulation domains.
 * For every simulation that has boundary conditions on this boundary
 * the properties/models are stored.
 */
class Boundary
{

  public:

    //! Constructor
    /*!
     * \param name a user defined name to identify this boundary
     */
    Boundary(const std::string& name);

    //! Destructor
    ~Boundary(void);

    //! Add a boundary property for simulation with ID \c simulator_id
    void add_boundary_properties(BoundaryProperties* properties,
        ID simulator_id);

    //! Get the boundary property for simulation with ID \c simulator_id
    BoundaryProperties* get_boundary_properties(ID simulator_id) const;

    //! Get the user defined name of this boundary
    const std::string& get_name(void) const;

    //! Initialize all models
    void init(void);


  private:

    //! A typedef for convenience
    typedef std::map<ID, BoundaryProperties*> PropertyMap;

    //! The models for the different simulations
    PropertyMap _models;

    //! The user defined name of this boundary
    std::string _name;

};

inline
BoundaryProperties*
Boundary::get_boundary_properties(ID simulator_id) const
{
  PropertyMap::const_iterator it(_models.find(simulator_id));

  if (it != _models.end())
    return it->second;
  else
    return NULL;
  
}

inline
Boundary::Boundary(const std::string& name)
  : _name(name)
{
}

inline
const std::string&
Boundary::get_name(void) const
{
  return _name;
}

#endif // _BOUNDARY_H_
