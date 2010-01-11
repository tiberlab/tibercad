#ifndef _OPTRECOMBINSPECTRUM_H_
#define _OPTRECOMBINSPECTRUM_H_






#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>


#include <fstream>
#include <iomanip>


#include <petsc_matrix.h>
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>

#include "EnvelopFunctionApprox.h"
#include "tensor.h"
#include "KspaceIntegration.h"

class OpticsKP;

//! K-space integration of  optical spectrum.

class OptRecombinSpectrum : public KspaceIntegration
{



 public:

  enum JobKind
  {
     RECOMBINATION  = 0, //<! spontaneous recombination spectrum
     ABSORPTION = 1 //<! optical absorbtion spectrum
  };

  struct options 
  {
    double Emin;//<! left boundary of spectrum [eV] 

    double Emax;//<! right boundary of spectrum [eV] 

    double dE;//<! spectrum mesh step [eV]

    double Gamma;//<! spectrum broadening [eV] 

    Tensor1 polariz; //<! light  polarization [eV] 


  };


  //!Constructor
  OptRecombinSpectrum(const ModelOptions& options);

  //! Destructor 
  virtual ~OptRecombinSpectrum();


  //!creates a new object 
  static  OptRecombinSpectrum* create(const ModelOptions& options);


 private:

  //!set options for the object
  void set_options( OptRecombinSpectrum::options& options   );


  //! pointers to   simulations which  solve Schroedinger equation 
  EnvelopFunctionApprox*  _quantum_model_initial_state;
  EnvelopFunctionApprox*  _quantum_model_final_state  ;

  //! pointer  to   simulation  which  calculates  optical matrix  element
  OpticsKP*  _optical_model;


  //!options 
  options opt;  


  //! Mesh for spectrum [eV];
  Mesh* _energy_mesh;

  //! spectrum type to calculate
  JobKind job;
  
 protected:


  virtual void 	do_init(void);

  virtual void 	parse_options (void);


  //!calculates objects k_point_density and eigen_energy
  virtual void calculate_for_k_point(const Point& k_point, 
				     std::map<const Elem*, double>& density, 
				     double& integrated_quantity);

  //! output of spectrum
  virtual void do_plot(void) ;


};


//---------------------------------------------------------


inline OptRecombinSpectrum*  OptRecombinSpectrum::create(const ModelOptions& options)
{
  return (new OptRecombinSpectrum(options) );
}

#endif
