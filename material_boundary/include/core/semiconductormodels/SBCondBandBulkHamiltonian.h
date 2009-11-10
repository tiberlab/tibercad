// $Id$


#ifndef _SBCONDBULKHAMILTONIAN_H_
#define _SBCONDBULKHAMILTONIAN_H_

#include "SBbulkHamiltonian.h"
#include "Semiconductor.h"


class Semiconductor;

//! A clas that builds a single band Hamiltonian for a conduction band of a crystal
class SBCondBandBulkHamiltonian: public SBbulkHamiltonian
{
 public:
  SBCondBandBulkHamiltonian(void);


 ~SBCondBandBulkHamiltonian(void);

  
 virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential) = 0;
  
  

 virtual void set_temperature(double Temperature);


 protected:

 virtual PhysicalModelInterface* create_new(void) const = 0 ;


 virtual void do_init(void);

 
 virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);


 //! a pointer to a semiconductor that contains parameters
 Semiconductor* semiconductor;

 //!calculate everything we need from semiconductor model
 virtual void calculate_for_init(void) = 0;

 private:


 

};
//----------------------------------------------------------------------------------------//

inline SBCondBandBulkHamiltonian::SBCondBandBulkHamiltonian()
  : semiconductor(NULL)
{
}



inline 
void SBCondBandBulkHamiltonian::set_temperature(double Temperature)
{
  semiconductor->set_temperature(Temperature);
}


#endif
