#ifndef _SBCONDBULKHAMILTONIAN_H_
#define _SBCONDBULKHAMILTONIAN_H_

#include "SBbulkHamiltonian.h"
#include "PhysicalModelInterface.h"
#include "Material.h"
//! A clas that builds a single band Hamiltonian for a conduction band of a crystal

class Semiconductor;

class SBCondBandBulkHamiltonian: public SBbulkHamiltonian
{
 public:
  SBCondBandBulkHamiltonian(void);


 ~SBCondBandBulkHamiltonian(void);

  
 virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential) = 0;
  
  

 protected:

 virtual PhysicalModelInterface* create_new(void) const = 0 ;


 virtual void do_init(void);

 
 virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);


 //! a pointer to a semiconductor that contains parameters
 Semiconductor* semiconductor;

 //!calculate everything we need from semiconductor model
 virtual void calculate_for_init(void) = 0;

 private:


 

};
//----------------------------------------------------------------------------------------//

inline SBCondBandBulkHamiltonian::SBCondBandBulkHamiltonian()
{
  semiconductor = NULL;
}

//-----------------------------------------------------------------------------------------//

inline SBCondBandBulkHamiltonian::~SBCondBandBulkHamiltonian()
{
  PhysicalModelInterface::destroy(semiconductor);
}

//-----------------------------------------------------------------------------------------//

inline void SBCondBandBulkHamiltonian::do_init()
{
  SBbulkHamiltonian::do_init();

  kp_bands.resize(1,0);
 
  kp_bands_map.insert(std::make_pair(0,0));


  const ModelOptions& opt =  get_options ();

  PhysicalModelInterface::destroy(semiconductor);

  semiconductor = Semiconductor::create( get_material() -> get_structure(), opt);

  semiconductor->init();
  
}

//------------------------------------------------------------------------------------------//

inline void SBCondBandBulkHamiltonian::calculate_VCA (const PhysicalModelInterface *comp_A, 
						      const PhysicalModelInterface *comp_B, double xa)
{
  const SBCondBandBulkHamiltonian* matA = dynamic_cast< const SBCondBandBulkHamiltonian*> (comp_A);

  const SBCondBandBulkHamiltonian* matB = dynamic_cast< const SBCondBandBulkHamiltonian*> (comp_B);

  semiconductor->build_alloy(matA->semiconductor, matB->semiconductor, xa);

  calculate_for_init();

}

#endif
