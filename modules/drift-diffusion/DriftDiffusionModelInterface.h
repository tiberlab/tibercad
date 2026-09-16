/*  
 * This file is part of the tiberCAD module driftdiffusion.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DriftDiffusionModelInterface.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_DRIFTDIFFUSIONMODELINTERFACE_H
#define TC_DRIFTDIFFUSIONMODELINTERFACE_H


#include "tibercad/base/TypeDefs.h"
#include "tibercad/physics/PhysicalModel.h"

#include <cassert>

class DriftDiffusionProperties;

//! Base class for the different models used in Drift-Diffusion calculations
/*!
 * This is the base class for all implementations of mobility models,
 * recombination models etc. which will be used in conjunction with
 * a semiconductor model derived from DriftDiffusionProperties
 */
class TC_DLEXPORT DriftDiffusionModelInterface : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~DriftDiffusionModelInterface(void) {};


    //! Get a reference to the DriftDiffusionProperties object
    /*!
     * \return a reference to the DriftDiffusionProperties object this
     * model belongs to
     */
    DriftDiffusionProperties& get_driftdiffusionproperties(void) const;


  protected:

    //! Empty constructor
    DriftDiffusionModelInterface(const ModelOptions& options);

    //! The standard reference temperature in eV
    /*!
     * We give this value here because we get the lattice temperature
     * from DriftDiffusionProperties in eV and many models need something
     * like T/T0
     */
    static const double T0;


  private:

    
};


inline
DriftDiffusionModelInterface::DriftDiffusionModelInterface(const ModelOptions& options)
  : PhysicalModel(options)
{
}






#endif // TC_DRIFTDIFFUSIONMODELINTERFACE_H
