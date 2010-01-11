// $Id$

#ifndef _FEMEIGENVALUEPROBLEM_H_
#define _FEMEIGENVALUEPROBLEM_H_

#include "EigenvalueProblem.h"

// Basic include files needed for the mesh functionality.
#include "mesh.h"
#include "linear_implicit_system.h"

#include "getpot.h"

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
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>


#include <fstream>
#include <iomanip>

//------------------------------------------------------------------------------

#include <complex>
#include <vector>
#include <petsc_matrix.h>
#include <algorithm>
#include <set>
#include "mesh_data.h"

class EquationSystems;

//! Base class to solve complex valued eigenvalue problems based on FEM
class FEMEigenvalueProblem : public  EigenvalueProblem
{
 public:

  //! To distinguish between H and S matrix
  enum MatrixName
  {
    H = 0,
    S = 1
  };

  //! data structure that contain information about if the dof is independent or not
  struct dof_new
  {
    bool independent;  //!< true if it is and independent dof
    unsigned int new_number;  //!< new number in the independent dofs list
  };

  enum Method
  {
    FEM = 0,
    BIM = 1
  };

  struct options
  {
    std::string solver;                 //!< solver type
   
    unsigned int max_iteration_number;  //!< maximum number of iterations for the eigenvalue solver
    
    double    eigen_solver_tolerance;   //!< tolerance for eigenvalue solver 
   
    bool solve_ev_problem_twice;        //!< if true, calculate the first eigenvalue only and then run again
   
    unsigned int number_of_eigenstates; //!< number of eigenstates to be calculated
    
    Method discretization_method;       //!< box integration or finite element
   
    bool Dirichlet_bc_everywhere;       //!< apply dirichlet boundary conditions at all boundaries

    std::string preconditioner;         //!< preconditioner name

    bool periodicity[3];                //!< periodic boundary conditions


    std::string st_ksp_type;            //!< Liner system solution method 
    

    std::string strategy;               //<! matlab (algorithm used in Matlab) or general (recommended by SLEPC)

    bool monitor;   //<! activates convergence monitor if true


    double spectrum_inversion_tolerance; //<! tolerance for spectrum inversion

    
  };


  FEMEigenvalueProblem(const ModelOptions& options);

  virtual ~FEMEigenvalueProblem() {};

 protected:

  virtual void do_init();

  virtual void do_solve() {};


  virtual void parse_options();



  //!calculates matricies H and S for the generalized problem Hx = gSx 
  virtual void calculate_Hamiltonian_and_S(void) {};

 


  //!dimension of the system
  short dim;

  //!pointer to the real part of the Hamiltonian
  SparseMatrix<Number>* Ham_real;

  //!pointer to the imaginary part of the Hamiltonian
  SparseMatrix<Number>* Ham_imag;

  //!pointer to the real part of S matrix 
  SparseMatrix<Number>* S_real;


  //!pointer to the real part of S matrix 
  SparseMatrix<Number>* S_imag;



  //!creates complex matrix inside EigenSolver class (TODO)
  /*!
    \param matrix  H or S
    \param m_real  real part of the matrix
    \param m_imag  imaginary part of the matrix. If NULL, matrix is considdred to be pure real
  */
  void create_complex_matrix(const MatrixName matrix, const SparseMatrix<Number>* m_real, 
  			     const SparseMatrix<Number>* m_imag = NULL ) {};


  //!number of independent dofs
  int number_of_new_dofs;


  //!total number of dofs 
  int number_of_all_dofs;



  //!vector: each element contains information about dof
  std::vector<FEMEigenvalueProblem::dof_new> new_dofs;


  //!creates new_dofs vector
  void make_new_dofs(void);


  //!pointer to mesh of the equation systems
  MeshBase* mesh;


   //!pointer the equation systems object
  EquationSystems* es;
  
  //!system that we add to the equation systems
  LinearImplicitSystem* system;


  //!diriclet DOFS
  std::set<unsigned int>  dirichlet_dofs;

  //!updates my_dof_constraints 
  void make_constraints(void);


  //! my copy of DofMap::_dof_constraints (I have to recalculate it because it is private)
  DofConstraints 	my_dof_constraints;

  
  //!creates dirichlet dofs
  void create_dirichlet_dofs(void);
  

  //!Apply Dirichlet boundary conditions to all boundary points!
  void apply_diriclet_bc_at_all_boundaries();
  

  //!put spectrum shift energy to be almost equal to the 1st eigenvalue
  virtual double get_new_spectrum_shift(void){};


 


  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spectrum shift 
  */
  void solve_eigen_value_problem(unsigned int ev_number, double spectrum_shift = 0.0 );


   //!read SLEPc solutions
  /*!
 
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want 
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) do something else
 
    \param number_of_ev number of eigen functions to read
  */
  virtual void read_SLEPC_solution(unsigned int number_of_ev) {};

  //!options of any eigensolver problem
  options solver_opt;



  //!passes H matrix to the eigensolver
  void copy_H_matrix_to_solver();



  //! Apply periodic boundary conditions
  void apply_periodic_bc();

  //! create list of nodes that lies at the periodic boundary
  void make_nodes_periodic();


  //!list of periodic nodes
  std :: vector< std :: vector <const Node*> >  nodes_periodic; //dim node list's: each contains list of nodes that periodic b.c
                                                                //must be applied to


  //!simulation domain boundary
  double min_coord[3];
  
  //!simulation domain boundary
  double max_coord[3];

  //! checks if element lies on boundary
  bool element_on_boundary(const Elem* element);
  



 private:

  //!size of the Hamiltonian matrix
  unsigned int _hamiltonian_size;


};


//---------------------------------------------------------------------------------//

inline bool FEMEigenvalueProblem::element_on_boundary(const Elem* element)
{
  bool result = false;

  
    
  unsigned int n_sides ; 

  if ( dim > 1 ) 
    n_sides = element->n_sides();
  else
    n_sides = element->n_nodes();


  for (short i = 0; i < n_sides; i++)
  {
    Elem* el1 = element->neighbor(i);

    if ( (el1 == NULL)  ) 
      result = true;
    else
      if (!( el1 -> active() ))
	result = true;
	  
    if (result) break;
	
      
  }

 
  return(result);
  
}


#endif // _FEMEIGENVALUEPROBLEM_H_
