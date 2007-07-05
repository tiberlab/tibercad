#ifndef _KSPACEINTEGRATION_H_
#define _KSPACEINTEGRATION_H_

#include "Kspace.h"

// Basic include files needed for the mesh functionality.

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


#include <fstream>
#include <iomanip>

#include "KspaceIntegration.h"
#include <map>


class KspaceIntegration : public Kspace
{
 public:

  //!options for charge density
  struct options
  {
    bool k_domain_user_input;       //!< if true user provides Brilluoin zone size. Otherwise, it is calculated by the program
    bool k_domain_refinement;       //!< if true, program will refine the Brilluoin zone adaptively
    bool uniform_refinement;        //!< if true, all the cells in the k-space are refined

    double refine_fraction;         //!< fraction of the elements to be refined
    double relative_accuracy;       //!< stop refinement if \f$ ||\rho_{i+1} - \rho_i||/||rho_i|| < \epsilon \f$
    unsigned int maximum_ref_level; //!< maximum level for k space refinement
   
    unsigned int degeneracy;        //!< degeneracy factor
  }; 


  
  KspaceIntegration();


  virtual ~KspaceIntegration();
  


 protected:

  
  //!calculates the quantity performing mesh refinement of k-space, if required.
  virtual void calculate_convergent_density(void);

  //!calculates density that is necessary for eack k-point and a number that will be used for refinement 
  virtual void calculate_for_k_point(const Point& k_point, 
				     std::map<const Elem*, double>& density, 
				     double& integrated_quantity) {};

  

  //!map from element in the k-grid to a real space density, which is a map between real space elements and density   
  std::map< const Elem*, std::map <const Elem*, double>  > kspace_local_density;



  //!map from Elem in the k-grid to an integrated quantity used for the refinement criterion
  std::map<const Elem*, double> kspace_integral;

   
  //!result after integration
  std::map<const Elem*, double> real_space_density;


  std::map<const Elem*, double> old_real_space_density; 



  virtual double estimate_error(void);

  
  //! equation system defined at k-space
  EquationSystems*  eq;
   
  //! system 
  LinearImplicitSystem* system;
   
   

   //!calculate density for a particular k-grid
   virtual void calculate_density();


  


   //!estimates error for mesh refinement KellyErrorEstimator is called
   void estimate_error_for_refinement(ErrorVector& error);


   //! size of the real_space_density vector
   unsigned int real_space_density_size;


   //!put data into opt
   virtual void parse_options(void);

   //!volumes of the elements
   std::map<const Elem*, double> volume;


  
   //!calculate  volumes of the elements
   void calculate_volumes(void);
   
   //!how many elementsin k-space has to be done
   unsigned int how_many_elements_to_do();
   

  
   options opt;

   
   virtual void do_solve(void);

   virtual void do_init(void);

 private:

   

};
#endif
