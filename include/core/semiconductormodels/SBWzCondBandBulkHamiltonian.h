#ifndef _SBWzCondBandBulkHamiltonian_h_
#define _SBWzCondBandBulkHamiltonian_h_


#include "SBbulkHamiltonian.h"
#include "WzDDsemiconductor.h"

//! A class to calculate single band Hamiltonian of wurtzite material
class SBWzCondBandBulkHamiltonian : public SBbulkHamiltonian
{

 public:
  //!Constructor
  /*!
    \param parameters parameters that describe a wurtzite crystal
  */
  SBWzCondBandBulkHamiltonian(WzDDsemiconductor::WzDDparameters& parameters);
  


  
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

 private:

  WzDDsemiconductor::WzDDparameters  par; 

};

#endif
