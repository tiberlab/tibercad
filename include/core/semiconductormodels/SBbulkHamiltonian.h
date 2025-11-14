#ifndef _SBBULKHAMILTONIAN_H_
#define _SBBULKHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"

//class DDsemiconductor;


class SBbulkHamiltonian : public EFAbulkHamiltonian
{

 public:

  //! destructor
  ~SBbulkHamiltonian(void);

 

 

  virtual void calculate_Hamiltonian_k_par(void);

 

  virtual void calculate_Hamiltonian_gen(void); 


  //! Applies ONLY potential. Strain in applied only in derived classes
  virtual void apply_strain_and_potential(Tensor2Gen& strain_crystal, double el_potential);

   





 

 protected:

  //! default constructor
  SBbulkHamiltonian(const ModelOptions& options);

  //!band edge
  double edge;

  //!\f$ \left(\frac{1}{m}\right)_{ij} \f$ tensor
  Tensor2Gen imass;

  //!The only matrix element of the Hamiltonian
  MatrixElement single_band_ham;

 
  //! a pointer to a semiconductor that contains parameters
  //DDsemiconductor* semiconductor;


  virtual PhysicalModel* create_new(void) const = 0;

  virtual void do_init(void) = 0;

  virtual void do_print_info(void);

 private:

  

  
 
};


#endif
