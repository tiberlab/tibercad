// $Id$

#ifndef _DYNAMICALMATRIX_H_
#define _DYNAMICALMATRIX_H_

#include "PhysicalModel.h"
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
   void get_dynamical_matrix(Tensor2Gen& D); 

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
   void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);

  Tensor2Gen _dynamical_matrix;
  

   PhononModel* _phonon_model;

};




inline
void 
DynamicalMatrix::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  // generates dynamical matrix in calculation system
  _dynamical_matrix = sym(RotMatrix * ( _dynamical_matrix * (RotMatrix.transpose())));


}

inline
void
DynamicalMatrix::get_dynamical_matrix(Tensor2Gen& dynamical_matrix) 
{
   this->re_init();
   dynamical_matrix = _dynamical_matrix;
}


#endif
