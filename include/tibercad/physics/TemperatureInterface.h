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
 * \file TemperatureInterface.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_TEMPERATUREINTERFACE_H
#define TC_TEMPERATUREINTERFACE_H

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/tiber_dll.h"


#include <vector>
#include <string>

class SimulationInterface;

namespace libMesh
{
class Elem;
class Point;
}


//! An interface to a temperature simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated temperature. If there is no simulation,
 * it will just return the temperature from the simulation options
 */
class TemperatureInterface
{

  public:

    //! Default constructor
    TemperatureInterface(void);


    //! Copy constructor
    TemperatureInterface(const TemperatureInterface& rhs);


    //! Assignment operator
    TemperatureInterface& operator=(const TemperatureInterface& rhs);


    //! Specify the temperature simulation to use
    /*!
     * Returns true if \c name refers to a valid temperature
     * simulation.
     */
    bool set_simulation(const std::string& name);


    //! Get the nodal temperatures
    void get_temperature(const libMesh::Elem* elem, std::vector<double>& temperatures);


    //! Get the temperatures in specified points
    void get_temperature(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
        std::vector<double>& temperatures, bool refcoord = false);


    //! Get the temperature in one point
    /*!
     * Set \c refcoord to true if p is given in reference coordinates
     */
    double get_temperature(const libMesh::Elem* elem, const libMesh::Point& p, bool refcoord = false);


    //! Tells if this interface has a simulation associated
    bool has_simulation(void) const;

   //! Tells if the  simulation is solved
    bool is_solved(void) const;

    //! Get the associated simulation
    SimulationInterface* get_simulation(void);



  private:

    //! The temperature simulation
    SimulationInterface* _simulation;


    //! The name of the temperature variable
    std::string _variable_name;


    //! The ID as returned from the simulation
    ID _id;


};


//
// inline members
//

inline
bool
TemperatureInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
TemperatureInterface::get_simulation(void)
{
  return _simulation;
}

#endif // TC_TEMPERATUREINTERFACE_H
