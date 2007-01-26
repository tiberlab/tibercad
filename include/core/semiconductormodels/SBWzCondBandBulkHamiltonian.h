#ifndef _SBWzCondBandBulkHamiltonian_h_
#define _SBWzCondBandBulkHamiltonian_h_


#include "SBbulkHamiltonian.h"
#include "WzSemiconductor.h"
#include "SBCondBandBulkHamiltonian.h"

//! A class to calculate single band Hamiltonian of wurtzite material
class SBWzCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:
  //!Constructor

  SBWzCondBandBulkHamiltonian(){};

  //!Destructor 
  ~SBWzCondBandBulkHamiltonian(){};


  
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

  
  static SBWzCondBandBulkHamiltonian* create();
  

 private:

 

  WzSemiconductor::WzDDparameters*  wz_par; 

 protected:


  virtual PhysicalModelInterface* create_new(void) const ;

 
  virtual void do_init(void);

};

inline PhysicalModelInterface* SBWzCondBandBulkHamiltonian::create_new() const
{
  return new SBWzCondBandBulkHamiltonian();
}

inline SBWzCondBandBulkHamiltonian* SBWzCondBandBulkHamiltonian::create()
{
  return new SBWzCondBandBulkHamiltonian();
}


#endif
