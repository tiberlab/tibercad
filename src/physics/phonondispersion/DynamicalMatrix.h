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
 * \file DynamicalMatrix.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_DYNAMICALMATRIX_H
#define TC_DYNAMICALMATRIX_H

#include "tibercad/physics/PhysicalModel.h"
#include "PhononModel.h"

class PhononModel;

//! Class to return the dynamical matrix of a lattice in k = Gamma
/*!

The lattice thermal conductivity must be in W/(cm K)

*/
class DynamicalMatrix : public PhysicalModel

{

public:
  
  //!Constructor 
   DynamicalMatrix(const ModelOptions& options);

   //!Destructor
  ~DynamicalMatrix(){};

   //!provides conductivity in simulation system W/(cm K)
   void get_dynamical_matrix(Tensor2& D); 

   virtual void re_init(void)=0;

   virtual void set_phonon_model(PhononModel* phonon_model) =0;

private:

 

protected:

  virtual void do_init (void)=0;

  virtual void read_database(void)=0;

  virtual void do_init_alloy (const PhysicalModel *comp_A,
      const PhysicalModel *comp_B, double xa); 

  virtual PhysicalModel* create_new (void) const =0;

  //!rotates dynamical matrix into the simulation system
   void rotate_to_calculation_system(const Tensor2& RotMatrix);

  Tensor2 _dynamical_matrix;
  

   PhononModel* _phonon_model;

};




inline
void 
DynamicalMatrix::rotate_to_calculation_system(const Tensor2& RotMatrix)
{

  // generates dynamical matrix in calculation system
  _dynamical_matrix = sym(RotMatrix * ( _dynamical_matrix * (RotMatrix.transpose())));


}

inline
void
DynamicalMatrix::get_dynamical_matrix(Tensor2& dynamical_matrix) 
{
   this->re_init();
   dynamical_matrix = _dynamical_matrix;
}


#endif
