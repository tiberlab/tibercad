#ifndef _SBUSERHAMILTONIAN_H_
#define _SBUSERHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "SBbulkHamiltonian.h"

//class DDsemiconductor;


class SBuserHamiltonian : public SBbulkHamiltonian
{

 public:
  //! default constructor
  SBuserHamiltonian(void);

  //! destructor
  ~SBuserHamiltonian(void);


  static  SBuserHamiltonian* create(void);

  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

 protected:



  virtual PhysicalModelInterface* create_new(void) const;

  virtual void copy_from (const PhysicalModelInterface *rhs);

  virtual void do_init(void);

 
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) {};

  


 private:

 

  
 
};

inline PhysicalModelInterface* create_new(void) 
{
  return new SBuserHamiltonian();
}





inline SBbulkHamiltonian* create()
{
  return new SBuserHamiltonian();
}
#endif
