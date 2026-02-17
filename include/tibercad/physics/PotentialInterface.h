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
 * \file PotentialInterface.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_POTENTIALINTERFACE_H
#define TC_POTENTIALINTERFACE_H

#include "tibercad/base/TypeDefs.h"


#include <set>
#include <vector>
#include <string>

class SimulationInterface;
class Atom;


namespace libMesh
{
  class Elem;
  class Point;
}


//! An interface to a Poisson simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated electrostatic potential.
 * The current implementation will return 0 whenever no data can be
 * obtained from the Potential simulation or if there is no
 * simulation at all.
 */
class PotentialInterface
{

public:

  //! Default constructor creates an empty interface that must be initialized by set_simulation
  PotentialInterface(void);

  //! Constructor that invokes initializations
  PotentialInterface(const std::string& name, const std::string& variable = "ElPotential");


  //! Specify the simulation to use 
  /*!
   * Returns true if \c name refers to a valid potential simulation
   * (examples of potential variables: ElPotential, eQFermi, hQFermi).
   */
  bool set_simulation(const std::string& name, const std::string& variable = "");


  //! Get the electrostatic potential in specified points
  void get_potential(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
      std::vector<double>& potentials, bool local_coord = false);


  //! Get the electrostatic potential in one point
  double get_potential(const libMesh::Elem* elem, const libMesh::Point& p,
      bool local_coord = false);


  //! Get the electrostatic potential at the position of an atom
  double get_potential(const Atom* atom);


  //! Tells if this interface has a simulation associated
  bool has_simulation(void) const;


  //! Get the associated simulation
  SimulationInterface* get_simulation(void);



private:

  //! The potential simulation
  SimulationInterface* _simulation;


  //! The ID as returned from the simulation
  ID _id;

};


//
// inline members
//

inline
bool
PotentialInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
PotentialInterface::get_simulation(void)
{
  return _simulation;
}


#endif // TC_POTENTIALINTERFACE_H
