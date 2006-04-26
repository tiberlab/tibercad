#ifndef _ENVELOPFUNCTIONAPPROX_H_ 
#define _ENVELOPFUNCTIONAPPROX_H_
//! A class that constructs Hamiltonian and S-matrix 

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "linear_implicit_system.h"
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

//------------------------------------------------------------------------------

#include <complex>
#include <vector>
#include <petsc_matrix.h>
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>
#include "mesh_data.h"

class EnvelopFunctionApprox
{
 public:
  //!control options


  //! data structure that contains options for effective mass
  struct options
  {
   
    unsigned int number_of_bands; //!< number of bands in EFA 

    std::string particle;   //!< particle name "el" or "hl"

    double  length_scale;   //!< mesh length scale [Bohr radius]

    double  max_energy;     //!< max energy value for Dirichlet b.c.

    bool periodicity[3];    //!< periodic boundary conditions

    std::string solver;  //!< solver type

    double eigen_solver_tolerance; //!< tolerance for eigenvalue solver [Ha]

    std::string output_type; //!< output type

  };





  //! data structure that contains eigenvalue and eigenvector
  struct eigen_propblem_solution
  {
    double eigen_energy; //!< eigen energy [eV]
    std::vector< std::complex<double>  > eigen_vector; //< eigen vector
  };


  //! data structure that contain information about if the dof is independent or not
  struct dof_new
  {
    bool independent;  //!< true if it is and independent dof
    unsigned int new_number;  //!< new number in the independent dofs list
  };


  //!constructor
  /*!
    \param opt  parameters  of the model
    \param mesh simulation domain mesh 
    \param mesh_data_input informations about materials and bondary conditions
  */
  EnvelopFunctionApprox(options& opt, Mesh& mesh, MeshData& mesh_data_input);


  //!destructor
  ~EnvelopFunctionApprox();


  //!computes Hamiltonian and S matrix
  void calculate_Hamiltonian_and_S(void); 



  //!sets material data
  void set_material_parameters(std::vector<EFAbulkHamiltonian*>&  bulkHamiltonian); 


  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param solution  eigenvectors and eigenvalues
  */
  void solve_eigen_value_problem(unsigned int ev_number );

  //!writes on disk the eigenfunction
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_eigen_functions(unsigned int state_number,  std::string& filename);


  //! assigned mesh_data_objects
  /*!
    \param  mesh_data_in mesh data that contains material information and Dirichlet conditions
   */
  void assign_mesh_data(MeshData& mesh_data_in);

  std::vector<eigen_propblem_solution> solution;


  //! define Dirichlet boundary condition
  /*!
    \param dirichlet_nodes_input nodes where \f$ \psivar({\bf r}) = 0 \f$
  */
  void define_diriclet_nodes(std::vector<unsigned int>&  dirichlet_nodes_input);


 private:

  options opt;
  
  MeshData*  meshdata;

  EquationSystems* es;

  //!diriclet nodes vector
  std::vector<unsigned int>  dirichlet_nodes;

  //!diriclet DOFS
  std::set<unsigned int>  dirichlet_dofs;
  
  //!my Jacobian because I calculate everything in atomic units
  double my_Jacobian; 

  
  //!dimension of the system
  short dim;

  //!pointer to the real part of the Hamiltonian
  SparseMatrix<Number>* Ham_real;

  //!pointer to the imaginary part of the Hamiltonian
  SparseMatrix<Number>* Ham_imag;

  //!pointer to the S matrix (it's real)
  SparseMatrix<Number>* S_real;

  //!vector that contains pointers to bulk Hamiltoninas
  std::vector<EFAbulkHamiltonian*>  bulkHamiltonian;

  //!bands names
  std::vector<std::string> psi_name;


  //!vector that contains material number
  std::vector<unsigned int> material_of_elem;


  //!create material_number vector
  void assemble_material_list(void);


  //!saves S matrix in PETSc format
  /*!
    \param file_name name of binary file for S matrix
  */
  void save_S_matrix(const std::string & file_name);

  //!saves H matrix in PETSc format
  /*!
    \param file_name name of binary file for S matrix
  */
  void save_H_matrix(const std::string & file_name);

  //!read SLEPc solutions
  void read_SLEPC_solution();

  //!vector: each element contains information about dof
  std::vector<EnvelopFunctionApprox::dof_new> new_dofs;

 
  //!swaps 4 byte variable for output
  inline void endian_swap(unsigned int& x)
    {
      x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
    }


  //!swaps 8 byte variable for output
  // __int64 for MSVC, "long long" for gcc
  inline void endian_swap(unsigned long long& x)
    {
      x = (x>>56) | 
        ((x<<40) & 0x00FF000000000000LL) |
        ((x<<24) & 0x0000FF0000000000LL) |
        ((x<<8)  & 0x000000FF00000000LL) |
        ((x>>8)  & 0x00000000FF000000LL) |
        ((x>>24) & 0x0000000000FF0000LL) |
        ((x>>40) & 0x000000000000FF00LL) |
        (x<<56);
    }

  
  //!creates dirichlet dofs
  void create_dirichlet_dofs(void);



  //!creates constraints
  void make_constraints(void);

  //!creates new_dofs vector
  void make_new_dofs(void)

};
#endif
