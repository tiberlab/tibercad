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
 * \file StrainInterface.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef _STRAININTERFACE_H_
#define _STRAININTERFACE_H_

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/tiber_dll.h"


#include <set>
#include <vector>
#include <string>

class SimulationInterface;
class Tensor2;

namespace libMesh
{
class Elem;
class Point;
}


//! An interface to a strain simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated strain. If there is no simulation,
 * it will just return 0 strain
 */
class StrainInterface
{

  public:

    //! Default constructor
    StrainInterface(void);


    //! Specify the strain simulation to use
    /*!
     * Returns true if name corresponds to an existing strain
     * simulation.
     */
    bool set_simulation(const std::string& name);

    
    //! Get the strain in calculation coordinate system
    void get_strain(const libMesh::Elem* elem, const libMesh::Point& point, Tensor2& strain);

    //! Get the strain in crystal coordinate system
    void get_crystal_strain(const libMesh::Elem* elem, const libMesh::Point& point, Tensor2& strain);

    //! Get the stress in calculation coordinate system
    void get_stress(const libMesh::Elem* elem, const libMesh::Point& point, Tensor2& stress);

    //! Get the stress in crystal coordinate system
    void get_crystal_stress(const libMesh::Elem* elem, const libMesh::Point& point, Tensor2& stress);


    //! Tells if this interface has a simulation associated
    bool has_simulation(void) const;


    //! Get the associated simulation
    SimulationInterface* get_simulation(void);


  private:

    //! The strain simulation
    SimulationInterface* _simulation;


    //! The ID for variable Strain
    ID _strain_id;

    //! The ID for variable StrainCrystal
    ID _strain_cryst_id;


    //! The ID for variable Stress
    ID _stress_id;

    //! The ID for variable StressCrystal
    ID _stress_cryst_id;

    

    //! Get a tensor
    void _get_data(const libMesh::Elem* elem, const libMesh::Point& point, Tensor2& data, ID id);
};


//
// inline members
//

inline
bool
StrainInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
StrainInterface::get_simulation(void)
{
  return _simulation;
}


#endif // _STRAININTERFACE_H_
