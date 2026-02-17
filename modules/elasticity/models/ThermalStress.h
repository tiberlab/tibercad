/*  
 * This file is part of the tiberCAD module elasticity.
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
 * \file ThermalStress.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_THERMALSTRESS_H
#define TC_THERMALSTRESS_H

#include "BodyForceModel.h"

#include "tibercad/physics/TemperatureInterface.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL ThermalStress : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ThermalStress(void);
  
  //! Creator function
  static ThermalStress* create(const ModelOptions& options);
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point);

  protected:

    //! Initialize
    virtual void do_init(void);


    //! Read expansion coefficients
    virtual void read_database(void);


  private:
  

    //! Constructor
    ThermalStress(const ModelOptions& options);

    //! Thermal expansion coefficients for the crystal directions
    libMesh::RealVectorValue _alpha;

    //! A reference temperature
    double _ref_temp;

    //! From where to get the temperature
    TemperatureInterface _temp;


  
};




inline
ThermalStress*
ThermalStress::create(const ModelOptions& options)
{ 
  return new  ThermalStress(options);
}




#endif // TC_THERMALSTRESS_H
