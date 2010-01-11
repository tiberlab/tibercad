#ifndef _SBZbCondBandBulkHamiltonian_h_
#define _SBZbCondBandBulkHamiltonian_h_


#include "ZbSemiconductor.h"
#include "SBCondBandBulkHamiltonian.h"
#include "PhysicalModelInterface.h"

//! A class to calculate single band Hamiltonian of zinc-blende material
class SBZbCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:

  //!Destructor
  ~SBZbCondBandBulkHamiltonian(){};
  
  //! Gamma or X or L
  std::string min_name;

  //! minima number 
  /*!
    O for Gamma; 
    0,1,2 for X;
    0,1,2,3 for L.
  */ 
  unsigned int min_number;

  
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

  static SBZbCondBandBulkHamiltonian* create(const ModelOptions& options);

 private:

  const ZbSemiconductor::ZbDDparameters* zb_par;

 protected:

  //!Constructor
  SBZbCondBandBulkHamiltonian(const ModelOptions& options)
    : SBCondBandBulkHamiltonian(options) {};

  virtual PhysicalModelInterface* create_new(void) const;


  virtual void calculate_for_init(void);


};

inline PhysicalModelInterface* SBZbCondBandBulkHamiltonian::create_new() const
{
  return new SBZbCondBandBulkHamiltonian(get_options());
}

inline SBZbCondBandBulkHamiltonian* SBZbCondBandBulkHamiltonian::create(const ModelOptions& options)
{
  return new SBZbCondBandBulkHamiltonian(options);
}

#endif
