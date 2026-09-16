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
 * \file Boundary.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_BOUNDARY_H
#define TC_BOUNDARY_H

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/ModelOptions.h"

#include <map>
#include <set>
#include <vector>
#include <string>

class PhysicalModel;
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

  private:

    //! A typedef for convenience
    typedef std::multimap<ID, PhysicalModel*> ModelMap;


  public:

    //! An iterator for the associated boundary models
    typedef ModelMap::iterator ModelIterator;

    //! A const iterator for the associated boundary models
    typedef ModelMap::const_iterator ConstModelIterator;


    //! Constructor
    /*!
     * \param name a user defined name to identify this boundary
     */
    Boundary(const std::string& names, const ModelOptions& options) TC_DLLOCAL;

    //! Destructor
    ~Boundary(void);

    //! Add a boundary property for simulation with ID \c simulator_id
    /*! \deprecated BoundaryProperties are deprecated */
    void add_boundary_properties(BoundaryProperties* properties,
        ID simulator_id) TC_DLLOCAL;

    //! Get the boundary property for simulation with ID \c simulator_id
    /*! \deprecated BoundaryProperties are deprecated */
    BoundaryProperties* get_boundary_properties(ID simulator_id) const;


    //! Add a boundary model and ID
    void add_model(ID id, PhysicalModel* model);


    //! Get the iterator to the first ID/model pair
    /*!
     * An ID might be occur several times, if boundaries with
     * different dimensions share the same ID.
     */
    ModelIterator models_begin(void);


    //! Get the past-the-end iterator for the model map
    ModelIterator models_end(void);


    //! Get the const iterator to the first ID/model pair
    /*!
     * An ID might be occur several times, if boundaries with
     * different dimensions share the same ID.
     */
    ConstModelIterator models_begin(void) const;


    //! Get the const past-the-end iterator for the model map
    ConstModelIterator models_end(void) const;


    //! Get the user defined name of this boundary
    const std::string& get_name(void) const;

    //! Initialize all models
    void init(void);

    //! Set the area factor
    void set_area_factor(double area_factor);

    //! Get the area factor
    double get_area_factor(void) const;

    //! Set the region IDs
    void set_region_ids(const std::vector<ID>& region_ids);

    //! Get the physical region IDs
    void get_region_ids(std::vector<ID>& ids) const;

    //! Get the region IDs as a set
    const std::set<ID>& get_region_ids(void) const;

    //! Check if it contains a region ID
    bool has_region_id(ID id) const;

    const ModelOptions& get_options(void) const;

  private:

    typedef std::map<ID, BoundaryProperties*> PropertyMap;

    //! The models for the different simulations
    ModelMap _models;

    //! Obsolete
    PropertyMap _oldmodels;

    //! The user defined name of this boundary
    std::string _name;

    //! The options
    ModelOptions _options;

    //! The area factor
    double _area_factor;

    //! The physical regions this boundary is touching
    std::set<ID> _region_ids;


};


//
// inline methods
//



inline
BoundaryProperties*
Boundary::get_boundary_properties(ID simulator_id) const
{
  PropertyMap::const_iterator it(_oldmodels.find(simulator_id));

  if (it != _oldmodels.end())
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
void
Boundary::add_model(ID id, PhysicalModel* model)
{
  _models.insert(std::make_pair(id, model));
}



inline
Boundary::ModelIterator
Boundary::models_begin(void)
{
  return _models.begin();
}



inline
Boundary::ModelIterator
Boundary::models_end(void)
{
  return _models.end();
}


inline
Boundary::ConstModelIterator
Boundary::models_begin(void) const
{
  return _models.begin();
}



inline
Boundary::ConstModelIterator
Boundary::models_end(void) const
{
  return _models.end();
}



inline
const std::set<ID>&
Boundary::get_region_ids(void) const
{
  return _region_ids;
}



inline
bool
Boundary::has_region_id(ID id) const
{
  return static_cast<bool>(_region_ids.count(id));
}

inline
const ModelOptions&
Boundary::get_options(void) const
{
  return _options;
}

#endif // TC_BOUNDARY_H
