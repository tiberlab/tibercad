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
 * \file ConstantBodyForce.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_CONSTANTBODYFORCE_H
#define TC_CONSTANTBODYFORCE_H

#include "BodyForceModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

//class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL ConstantBodyForce : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ConstantBodyForce(void) {};
  
  //! Creator function
  static ConstantBodyForce* create(const ModelOptions& options);
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point){};

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModel* comp_A,
    //         const PhysicalModel* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);



  private:
  

  //! Constructor
     ConstantBodyForce(const ModelOptions& options);
  
};




inline
ConstantBodyForce*
ConstantBodyForce::create(const ModelOptions& options)
{ 
  return new  ConstantBodyForce(options);
}




#endif // TC_GRAYMODEL_H
