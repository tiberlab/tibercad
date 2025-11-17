#ifndef _SBBULKHAMILTONIAN_H_
#define _SBBULKHAMILTONIAN_H_

#include "EFAbulkHamiltonian.h"
#include "Tensor2.h"

#include <complex>
#include <vector>



//! A class that builds single band Hamiltonian
class SBbulkHamiltonian : public EFAbulkHamiltonian
{

 public:

  //! destructor
  ~SBbulkHamiltonian(void);

 

 

  virtual void calculate_Hamiltonian_k_par(void);

 

  virtual void calculate_Hamiltonian_gen(void); 


  //! Applies ONLY potential. Strain in applied only in derived classes
  virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential);

   





 

 protected:

  //! default constructor
  SBbulkHamiltonian(const ModelOptions& options);

  //!band edge
  double edge;

  //!\f$ \left(\frac{1}{m}\right)_{ij} \f$ tensor
  Tensor2 imass;

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
