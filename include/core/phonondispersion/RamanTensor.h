// $Id$

#ifndef _RAMANTENSOR_H_
#define _RAMANTENSOR_H_

#include "PhysicalModelInterface.h"
#include "PhononModel.h"

class PhononModel;


//! Class to return the dynamical matrix of a lattice in k = Gamma
class RamanTensor : public PhysicalModelInterface
{

public:
  
  //!Constructor 
   RamanTensor();

   //!Destructor
  ~RamanTensor(){};

   void get_raman_tensor(std::vector< Tensor2Sym> & RT); 

   virtual void re_init(void)=0;

   virtual void set_phonon_model(PhononModel* phonon_model) =0;

private:

 

protected:

  virtual void do_init (void)=0;

  virtual void read_database(void)=0;

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A,
      const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const =0;

  //!rotates dynamical matrix into the simulation system
   void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);

  std::vector< Tensor2Sym >  _raman_tensor;
  
  PhononModel* _phonon_model;

};



inline
void 
RamanTensor::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  // generates dynamical matrix in calculation system
   for (unsigned int n= 0; n<3; n++)
{  _raman_tensor[n] = sym(RotMatrix * ( _raman_tensor[n] * (RotMatrix.transpose())));}


}

inline
void
RamanTensor::get_raman_tensor(std::vector<Tensor2Sym>& raman_tensor) 
{
   //this->re_init();
   raman_tensor = _raman_tensor;
}


#endif
