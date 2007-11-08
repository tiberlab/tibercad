#ifndef _QUANTUMDISPERSION_H_ 
#define _QUANTUMDISPERSION_H_

//! A class that calculates quantum density, performing integration of the density in k-space

// Basic include files needed for the mesh functionality.


#include "EnvelopFunctionApprox.h"
#include "Kspace.h"

//! This class is to calculate dispersion of a quantized state
class QuantumDispersion : public Kspace
{




 public:
  
  struct options
  {
    unsigned int min_eigenvalue_number;
    unsigned int max_eigenvalue_number;
    bool bulk_calculation;
  };
  

  

  //! Consructor
  QuantumDispersion();

  //! Desructor
  ~QuantumDispersion();


   
  //!creates a new object 
  static  QuantumDispersion* create();

 private:

  options opt;
  
  //! name of the simulation that solves Schroedinger equation 
  EnvelopFunctionApprox*  quantum_model;
 

  //! copy-constructor should not be used
  QuantumDispersion(const QuantumDispersion& t ) {}; 



  //!spectra of eigenvalues
  /*!
    eigen_energy[kspace_node_id][eigen_value_number];
  */
  std::vector<   std::vector<double>  >  eigen_energy;


  //!calculate eigen_energy
  void calculate_eigen_energy();

 protected:

  

 
   virtual void 	do_init(void);

   virtual void 	do_solve (void);

   virtual void 	parse_options (void);

   //!is reimplemented here to be able to do output over k-space 
   virtual void 	do_plot (void);
};


//---------------------------------------------------------

inline QuantumDispersion*  QuantumDispersion::create()
{
  return (new QuantumDispersion );
}

#endif
