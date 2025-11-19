#ifndef _SBWzCondBandBulkHamiltonian_h_
#define _SBWzCondBandBulkHamiltonian_h_


#include "tibercad/physics/semiconductormodels/SBbulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/WzSemiconductor.h"
#include "tibercad/physics/semiconductormodels/SBCondBandBulkHamiltonian.h"

//! A class to calculate single band Hamiltonian of wurtzite material
class SBWzCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:

  //!Destructor 
  ~SBWzCondBandBulkHamiltonian(){};


  
  virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential);

  
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
