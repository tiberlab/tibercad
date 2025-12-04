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
 * \file BuiltInStrain.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef _BUILTINSTRAIN_H_
#define _BUILTINSTRAIN_H_

#include "BodyForceModel.h"

#include "tensor_value.h"




//! The base class for Poisson boundary conditions
class TBDLLOCAL BuiltInStrain : public BodyForceModel
{

  public:

    //! Destructor
    ~BuiltInStrain(void) {};

    //! Creator function
    static BuiltInStrain* create(const ModelOptions& options);


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
    BuiltInStrain(const ModelOptions& options);

   
};



inline
BuiltInStrain::BuiltInStrain(const ModelOptions& options) :
  BodyForceModel(options)
{
}



inline
BuiltInStrain*
BuiltInStrain::create(const ModelOptions& options)
{
  return new BuiltInStrain(options);
}




#endif // _POISSONDIRICHLET_H_
