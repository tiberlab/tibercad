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
#include "SimulationInterface.h"
#include "tensor.h"

//! This class is a nextnano-like model of quantum density. The approach is correct, in principle, only for equilibrium.
class QuantumDensity : public SimulationInterface
{




 public:
  
  enum JobKind
  {
     DENSITY = 0,
     DISPERSION = 1
  };
  

  //!options for charge density
  struct options
  {
    bool k_domain_user_input;       //!< if true user provides Brilluoin zone size. Otherwise, it is calculated by the program
    bool k_domain_refinement;       //!< if true, program will refine the Brilluoin zone adaptively
    bool uniform_refinement;        //!< if true, all the cells in the k-space are refined

    double refine_fraction;         //!< fraction of the elements to be refined
    double relative_accuracy;       //!< stop refinement if \f$ ||\rho_{i+1} - \rho_i||/||rho_i|| < \epsilon \f$
    unsigned int maximum_ref_level; //!< maximum level for k space refinement

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



  //!returns reference to kmesh object
  const Mesh& get_k_mesh(void) const; 



  //!returns \f$ \rho({\bf k} ) = \int \rho{\bf{ k, r}}  d{\bf r} \f$
   std::vector<double>  get_density_in_k_space(void)  const;

   //!creates a new object 
   static  QuantumDensity* create();

 private:

     //!defines 1D  Brilluoin zone \f$ k \in [-{\bf K}/2; {\bf K}/2) \f$
  /*!
    \param k_vector Basis k-vector  \f$ \bf K \f$ [atom. units]
    \param n - initial number of nodes
  */
  void define_k_space(Tensor1 k_vector, unsigned int n);

 //!defines 2D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2) \f$
  /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param n - initial number of nodes in direction 1
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param m - initial number of nodes in direction 2
  */
  
  void define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m);

  //!defines 2D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2)   \otimes [-{\bf K}_3 / 2; {\bf K}_3 / 2) \f$
  /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param n - initial number of nodes in direction 1
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param m - initial number of nodes in direction 2
    \param k_vector3  Basis k-vector  \f$ {\bf K}_3 \f$ [atom. units]
    \param k - initial number of nodes in direction 3
  */
  
  void define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, unsigned int m, Tensor1 k_vector3, unsigned int k);

  

  


  //!set options for the object
  void set_options( QuantumDensity::options& options   );


 

  //!calculates density performing mesh refinement of k-space, if required.
  void calculate_convergent_density(void);



   //! name of the simulation that solves Schroedinger equation 
   EnvelopFunctionApprox*  quantum_model;


   //!quantum density for elements 
   std::map<const Elem*, double> real_space_density;


   //! size of the real_space_density vector
   unsigned int real_space_density_size;


   //!Brilluoin zone
   Mesh* kmesh;
  
   //!Boundaries of the Brilluoin zone [atomic units]
   double kmin[3], kmax[3];


   //!Dimension of the Brilluoin zone
   short  k_dim;


   //options 
   options opt;  

   //number of nodes in k-domain
   std::vector<unsigned int>  num_nodes;


   //!build k space grid
   void build_k_grid();

   //!map from node in the k-grid to a real space density, which is a map between real space elements and density   
   std::map< const Node*, std::map <const Elem*, double>  > k_point_density;

   //!map from node in the k-grid to a total charge
   std::map< const Node*, double > k_point_charge;


   //!spectra of eigenvalues
   /*!
     eigen_energy[i][j] = E_i(k_j);
   */
   std::vector< std::vector <double> > eigen_energy;

   

   //! Rotate mesh
   /*!
     \param mesh  pointer to the mesh
     \RotMatrix transformation matrix (not necessaryly rotation matrix) 
   */
   void rotate_mesh(Mesh* mesh, Tensor2Gen& RotMatrix);


   //! matrix that rotates mesh
   Tensor2Gen transform_matrix;

   //! equation system defined at k-space
   EquationSystems*  eq;
   
   //! system 
   LinearImplicitSystem* system;
   
   

    //!calculate density for a particular k-grid
   void calculate_density();


   //!put charge of each k-point into system solution
   void prepare_system_solution();


   //!calculates objects k_point_density and eigen_energy
   void calculate_at_each_k_point();
   

 protected:

   //!in this class  outputs particle density in a real space
   /*!
     Significant variable for this class is "quantum_density"
   */
   virtual void build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);

 
   virtual void 	do_init(void);

   virtual void 	do_solve (void);

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
