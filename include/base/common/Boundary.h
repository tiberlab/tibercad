// $Id$

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "TypeDefs.h"

#include <map>
#include <set>
#include <string>

class BoundaryProperties;
class SimulationEnvironment;


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
     * \param environment the environment
     * \param region_ids the boundary region IDs
     */
    Boundary(const std::string& name, SimulationEnvironment* environment,
        std::set<ID> region_ids);

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

    //! Set the area factor
    void set_area_factor(double area_factor);

    //! Get the area factor
    double get_area_factor(void) const;

    //! Get the simulation environment
    SimulationEnvironment* get_environment(void);

    //! Get the physical region IDs
    const std::set<ID>& get_region_ids(void) const;



  private:

    //! A typedef for convenience
    typedef std::map<ID, BoundaryProperties*> PropertyMap;

    //! The models for the different simulations
    PropertyMap _models;

    //! The user defined name of this boundary
    std::string _name;

    //! The area factor
    double _area_factor;

    //! The environment this boundary is belonging to
    SimulationEnvironment* _env;

    //! The physical regions this boundary is touching
    std::set<ID> _region_ids;

    //! Find the physical region IDs
    void find_region_ids(void);

};


//
// inline methods
// 



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
const std::string&
Boundary::get_name(void) const
{
  return _name;
}


inline
void
Boundary::set_area_factor(double area_factor)
{
  _area_factor = area_factor;
}


inline
double
Boundary::get_area_factor(void) const
{
  return _area_factor;
}


inline
const std::set<ID>&
Boundary::get_region_ids(void) const
{
  return _region_ids;
}


inline
SimulationEnvironment*
Boundary::get_environment(void)
{
  return _env;
}



#endif // _BOUNDARY_H_
