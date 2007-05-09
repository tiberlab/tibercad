#ifndef _QUANTUMDENSITYCALCULATION_H_ 
#define _QUANTUMDENSITYCALCULATION_H_

//! A class that calculates quantum density, performing integration of the density in k-space

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "equation_systems.h"

#include "getpot.h"
// For mesh refinement


#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"

// Define the Finite Element object.
#include "fe.h"
#include "elem.h"
// Define Gauss quadrature rules.
#include "quadrature_gauss.h" 

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

#include "fe_interface.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>


#include <fstream>
#include <iomanip>



#include <complex>
#include <vector>

#include <petsc_matrix.h>
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>
#include "mesh_data.h"
#include "DriftDiffusion.h"
#include "EnvelopFunctionApprox.h"
#include "tensor.h"
#include "KspaceIntegration.h"

//! This class is a nextnano-like model of quantum density. The approach is correct, in principle, only for equilibrium.
class QuantumDensity : public KspaceIntegration
{




 public:
  
  
  

  //!options for charge density
  struct options
  {
    

    double Temperature;             //!< temperature [K]
    unsigned int degeneracy;        //!< degeneracy factor to mutiply the charge density 

    unsigned int intial_eigenstates_number;  //!< number of required eigenstates for the first call of Schoedinger solver

  
    bool log_output; 
   
   

  };



  //! Consructor
  QuantumDensity();

  //! Desructor
  ~QuantumDensity();
  


  //!returns charge density for quadratur points 
  /*!
    We assume that the charge density is constant of the the element of a mesh used for a quantum calculation
    \param element pointer to the element
    \param quad_points quadratur points in local coordinates
    \param density resulting density
  */
  void get_particle_density(const Elem* element, const std::vector<double>& quad_points, std::vector<double> density);



  //!returns \f$ \rho({\bf k} ) = \int \rho{\bf{ k, r}}  d{\bf r} \f$
  std::vector<double>  get_density_in_k_space(void)  const;

  //!creates a new object 
  static  QuantumDensity* create();

 private:



  //!set options for the object
  void set_options( QuantumDensity::options& options   );




   //! name of the simulation that solves Schroedinger equation 
   EnvelopFunctionApprox*  quantum_model;


  


 
   //options 
   options opt;  



   //!calculates objects k_point_density and eigen_energy
   virtual void calculate_at_each_k_point();


  

 protected:

   //!in this class  outputs particle density in a real space
   /*!
     Significant variable for this class is "quantum_density"
   */
   virtual void build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);

 
   virtual void 	do_init(void);

   virtual void 	parse_options (void);

   //!is reimplemented here to be able to do output over k-space 
   virtual void 	do_plot (void);
};


//---------------------------------------------------------

inline QuantumDensity*  QuantumDensity::create()
{
  return (new QuantumDensity );
}

#endif
