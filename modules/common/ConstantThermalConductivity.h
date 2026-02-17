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
 * \file ConstantThermalConductivity.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_CONSTANTTHERMALCONDUCTIVITY_H
#define TC_CONSTANTTHERMALCONDUCTIVITY_H

#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/physics/misc/ThermalConductivityModel.h"





//! The base class for Poisson boundary conditions
class TBDLLOCAL ConstantThermalConductivity : public ThermalConductivityModel
{
  
public:
  
  //! Destructor
  virtual ~ConstantThermalConductivity(void) {};
  
  //! Creator function
  static ConstantThermalConductivity* create(const ModelOptions& options);
  
  virtual void calculate(const Elem* elem, const Point& point, double temperature);
  
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
  virtual void read_database(void);
  
  
  /* We do not use this here: */
  //  virtual void read_interface_database(void);
  
  
  //! Create a new object of the same type
  virtual PhysicalModel* create_new(void) const;
  
  
private:
  
  //! The thermal conductivity at reference temperature
  libMesh::RealGradient _kappa;

  //! The temperature coefficient
  libMesh::RealGradient _temp_coeff;

  //! The reference temperature (default 300K)
  double _ref_temp;

  //! Constructor
  ConstantThermalConductivity(const ModelOptions& options);
  
};


inline
PhysicalModel*
ConstantThermalConductivity::create_new(void) const
{
  return new   ConstantThermalConductivity(get_options());
}

inline
ConstantThermalConductivity*
ConstantThermalConductivity::create(const ModelOptions& options)
{
  return new  ConstantThermalConductivity(options);
}



#endif // TC_GRAYMODEL_H

  
