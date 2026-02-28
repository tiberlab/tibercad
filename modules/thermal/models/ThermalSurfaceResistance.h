/*  
 * This file is part of the tiberCAD module thermal.
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
 * \file ThermalSurfaceResistance.h
 * \brief tiberCAD thermal module header.
 *
 * \note This file is part of module thermal.
 */


#ifndef TC_THERMALSURFACERESISTANCE_H
#define TC_THERMALSURFACERESISTANCE_H

#include "ThermalBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"


namespace libMesh
{
  class Elem;
}



//! The base class for Poisson boundary conditions
class TC_DLLOCAL ThermalSurfaceResistance : public ThermalBoundaryModel
{

  public:

    //! Destructor
    ~ThermalSurfaceResistance(void) {};

    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) override;


  protected:

    //! Constructor
    explicit ThermalSurfaceResistance(const ModelOptions& options);

    //! Initialize
    virtual void do_init(void);


  private:

    double _temperature;
    double _resistance;
};



inline
ThermalSurfaceResistance::ThermalSurfaceResistance(const ModelOptions& options) :
  ThermalBoundaryModel(options),
  _temperature(0),
  _resistance(0)
{
}





#endif // TC_POISSONDIRICHLET_H
