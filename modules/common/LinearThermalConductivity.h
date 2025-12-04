/*  
 * This file is part of the tiberCAD module common.
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
 * \file LinearThermalConductivity.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _LINEARTHERMALCONDUCTIVITY_H_
#define _LINEARTHERMALCONDUCTIVITY_H_

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"

#include "tibercad/base/tiber_dll.h"

#include "elem.h"
#include "tibercad/physics/misc/ThermalConductivityModel.h"






//! The base class for Poisson boundary conditions
class TBDLLOCAL LinearThermalConductivity : public ThermalConductivityModel
{
  
public:
  
  //! Destructor
  virtual ~LinearThermalConductivity(void) {};
  
  //! Creator function
  static LinearThermalConductivity* create(const ModelOptions& options);
  
 virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double temperature);
  
protected:
  
  //! Initialize
  virtual void do_init(void);
  
  /* In some cases it might be useful to reimplement this: */
  // virtual void do_init_interface(const PhysicalModel* comp_A,
  //         const PhysicalModel* comp_B);
  
  virtual void do_init_alloy (const PhysicalModel *comp_A,
			      const PhysicalModel *comp_B, double xa){};
  
  virtual void  read_database_alloy(void){};
  /* This is not used here: */
  virtual void read_database(void){};
  
  
  /* We do not use this here: */
  //  virtual void read_interface_database(void);
  
  
  
private:
   
   double kx0;
   double kz0;
   double mx;
   double mz;
   double z0;
  
  //! Constructor
  LinearThermalConductivity(const ModelOptions& options);
  
};




inline
LinearThermalConductivity*
LinearThermalConductivity::create(const ModelOptions& options)
{
  return new  LinearThermalConductivity(options);
}




#endif // _GRAYMODEL_H_
