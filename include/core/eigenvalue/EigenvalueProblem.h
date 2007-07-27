#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_
//---------------------------------------------------------
//!Calculates optical modes 
// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "equation_systems.h"

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
#include "dense_submatrix.h"
#include "dense_subvector.h"
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
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>
#include "mesh_data.h"

//------------------------------------------------------------------------------
#include "SimulationInterface.h"

//! Abstract class to solve complex valued eigenvalue problem

class EigenvalueProblem: public  SimulationInterface
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

  EigenvalueProblem(void) {};

  virtual ~EigenvalueProblem() {};

 protected:

  virtual void do_init() {};

  virtual void do_solve() {};


  virtual void parse_options() {};


  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException) {};
    
   
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException) {};


  //!calculates matricies H and S for the generalized problem Hx = gSx 
  virtual void calculate_Hamiltonian_and_S(void) {};

  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spetrum shift [work units]
  */
  virtual void 	solve_eigen_value_problem (unsigned int ev_number, double spectrum_shift=0.0) {};


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



  //!creates complex matrix inside EigenSolver class
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
  std::vector<EigenvalueProblem::dof_new> new_dofs;

};

#endif 
