#ifndef _SBWzCondBandBulkHamiltonian_h_
#define _SBWzCondBandBulkHamiltonian_h_


#include "SBbulkHamiltonian.h"
#include "WzSemiconductor.h"
#include "SBCondBandBulkHamiltonian.h"

//! A class to calculate single band Hamiltonian of wurtzite material
class SBWzCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:

  //!Destructor 
  ~SBWzCondBandBulkHamiltonian(){};


  
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

  
  static SBWzCondBandBulkHamiltonian* create(const ModelOptions& options);
  

 private:

 

  const WzSemiconductor::WzDDparameters*  wz_par; 

 protected:

  //!Constructor
  SBWzCondBandBulkHamiltonian(const ModelOptions& options)
    : SBCondBandBulkHamiltonian(options) {};

  virtual PhysicalModel* create_new(void) const ;

 
  virtual void calculate_for_init(void);


};

inline PhysicalModel* SBWzCondBandBulkHamiltonian::create_new() const
{
  return new SBWzCondBandBulkHamiltonian(get_options());
}

inline SBWzCondBandBulkHamiltonian* SBWzCondBandBulkHamiltonian::create(const ModelOptions& options)
{
  return new SBWzCondBandBulkHamiltonian(options);
}


#endif
