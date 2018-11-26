#ifndef _SBUSERHAMILTONIAN_H_
#define _SBUSERHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "SBbulkHamiltonian.h"
#include "PhysicalModelInterface.h"
//class DDsemiconductor;


class SBuserHamiltonian : public SBbulkHamiltonian
{

 public:

  //! destructor
  ~SBuserHamiltonian(void);


  static  SBuserHamiltonian* create(const ModelOptions& options);

  

 protected:

  //! default constructor
  SBuserHamiltonian(const ModelOptions& options);



  virtual PhysicalModelInterface* create_new(void) const;

 

  virtual void do_init(void);

 
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) {};

  


 private:

 

  
 
};

inline PhysicalModelInterface* SBuserHamiltonian::create_new() const
{
  return new SBuserHamiltonian(get_options());
}





inline SBuserHamiltonian* SBuserHamiltonian::create(const ModelOptions& options)
{
  return new SBuserHamiltonian(options);
}
#endif
