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
 * \file BoundaryProperties.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_BOUNDARYPROPERTIES_H
#define TC_BOUNDARYPROPERTIES_H


#include "tibercad/module/TiberModelObject.h"
#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/TypeDefs.h"

class Boundary;

//! Boundary properties for a certain type of simulation.
/*!
 *  This is the base class for boundary properties. Every solver
 *  module should implement its own class derived from this one to hold
 *  the properties needed for calculations.
 */
class BoundaryProperties : public TiberModelObject
{

  public:

    //! The empty destructor
    /*!
     * should be implemented in the derived classes if needed
     */
    virtual ~BoundaryProperties(void);

    //! Initialize this boundary
    void init(void);

    //! Set the boundary this model is associated with
    void set_boundary(Boundary* boundary);

    //! Set the simulation ID
    void set_simulation_id(ID simulation_id);


  protected:

    //! The empty constructor
    /*!
     * \c BoundaryType should not be instantiated directly
     */
    BoundaryProperties(const ModelOptions& options);


    //! Initialize the model
    /*!
     * This method should set all model options and call
     * \c init() of any associated model
     *
     * This method should be implemented in derived classes.
     */
    virtual void do_init(void) {};

    //! Get the boundary this model is associated with
    Boundary* get_boundary(void);

    //! Get the simulation ID
    ID get_simulation_id(void) const;


  private:

    //! Disable copy constructor
    BoundaryProperties(const BoundaryProperties&);

    //! Disable assignement operator
    BoundaryProperties& operator=(const BoundaryProperties&);

    //! The boundary this model is associated with
    Boundary* _boundary;

    //! The ID of the simulation this model is associated with
    ID _simulation_id;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
BoundaryProperties::BoundaryProperties(const ModelOptions& options)
 : TiberModelObject(options)
{
}


inline
void
BoundaryProperties::init(void)
{
  do_init();
}



inline
void
BoundaryProperties::set_boundary(Boundary* boundary)
{
  _boundary = boundary;
}



inline
Boundary*
BoundaryProperties::get_boundary(void)
{
  return _boundary;
}


inline
void
BoundaryProperties::set_simulation_id(ID simulation_id)
{
  _simulation_id = simulation_id;
}


inline
ID
BoundaryProperties::get_simulation_id(void) const
{
  return _simulation_id;
}


#endif // TC_BOUNDARYPROPERTIES_H
