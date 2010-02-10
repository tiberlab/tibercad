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
#include "SimulationOptions.h"
#include "EnvelopFunctionApprox.h"
#include "tensor.h"
#include "KspaceIntegration.h"

//! This class is a nextnano-like model of quantum density. The approach is correct, in principle, only for equilibrium.
class QuantumDensity : public KspaceIntegration
{




 public:
  
  enum Variables
  {
    DENSITY = 0
  };
  

  //!options for charge density
  struct options
  {
    
    bool analitic; //!estimate density analitically rather than calculte numerically
    


   
    unsigned int degeneracy;        //!< degeneracy factor to mutiply the charge density 

    unsigned int intial_eigenstates_number;  //!< number of required eigenstates for the first call of Schoedinger solver

  
    bool log_output; //!<perform some screen output for debugging


    bool bulk_calculation; //!< bulk charge density 
   
   

  };



  //! Consructor
  QuantumDensity(const ModelOptions& options);

  //! Desructor
  ~QuantumDensity();
  


  //!returns charge density for quadratur points 
  /*!
    We assume that the charge density is constant of the the element of a mesh used for a quantum calculation
    \param element pointer to the element
    \param quad_points quadratur points in local coordinates
    \param density resulting density
  */
  void get_particle_density(const Elem* element, const std::vector<Point>& quad_points, std::vector<double>& density);



 
 

  //!creates a new object 
  static  QuantumDensity* create(const ModelOptions& options);



  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;
     
     
 
  virtual void get_solution_secure(const Elem* elem,
         const std::set<ID>& ids, std::vector<std::map<ID, double> >& values); 

 
  virtual void get_solution_secure(const Elem* elem,
				   const std::vector<Point>& p, const std::set<ID>& ids,
				   std::vector<std::map<ID, double> >& values);
  




  //!returns particle charge
  inline double get_particle_charge() const;

 private:



  //!set options for the object
  void set_options( QuantumDensity::options& options   );




   //! name of the simulation that solves Schroedinger equation 
   EnvelopFunctionApprox*  quantum_model;


  


 
   //options 
   options opt;  



   //!analitic (parabolic)  density calculation
   void estimate_analitic_density(void) ;
  

   //!k vector for effective mass calculation
   std::vector<double> k_vector1;

   //!k vector for effective mass calculation
   std::vector<double> k_vector2;

   //!vector of solution necessary for self-consistent solver
   AutoPtr< NumericVector<double> >  _solution_vector;


 protected:

   //!in this class  outputs particle density in a real space
   /*!
     Significant variable for this class is "quantum_density"
   */
   virtual void build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);

 
   virtual void 	do_init(void);

   virtual void 	parse_options (void);



   //!calculates objects k_point_density and eigen_energy
   virtual void calculate_for_k_point(const Point& k_point, 
				     std::map<const Elem*, double>& density, 
				     double& integrated_quantity);



   
   virtual void do_solve();


   virtual NumericVector< double > & 	do_get_solution_vector (void);

   //!here it copies new_solution both into _solution_vector and KspaceIntegration::real_space_density
   virtual void do_set_solution_vector(const NumericVector<double> & new_solution);


};


//---------------------------------------------------------

inline QuantumDensity*  QuantumDensity::create(const ModelOptions& options)
{
  return (new QuantumDensity(options) );
}

//---------------------------------------------------------

inline double QuantumDensity::get_particle_charge() const
{
  return( quantum_model->get_particle_charge() );
}

//---------------------------------------------------------
inline NumericVector< double > &  QuantumDensity::do_get_solution_vector (void)
{

 

  return (*_solution_vector);

}
//--------------------------------------------------------
#endif
