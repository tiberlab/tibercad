#ifndef _ENVELOPFUNCTIONAPPROX_H_ 
#define _ENVELOPFUNCTIONAPPROX_H_
//! A class that constructs Hamiltonian and S-matrix 

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

//------------------------------------------------------------------------------

#include <complex>
#include <vector>
#include <petsc_matrix.h>
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>
#include "mesh_data.h"
#include "macrostrain.h"
#include "DriftDiffusion.h"
class EnvelopFunctionApprox
{
 public:
  //!control options


  //! data structure that contains options for effective mass
  struct options
  {
    
    unsigned int number_of_bands; //!< number of bands in EFA 

    std::map<short, short> kp_bands; //!< map between band numbers: from 8 band scheem to any band scheem


    std::string particle;   //!< particle name "el" or "hl"

    double  length_scale;   //!< mesh length scale [Bohr radius]

    double  max_energy;     //!< max energy value for Dirichlet b.c.

    bool periodicity[3];    //!< periodic boundary conditions

    std::string solver;  //!< solver type

    std::string mpi_command_line; //!< something like:  mpirun -np 5 -machinefile machines

    std::string solver_command_line; //!< additional line to pass parameters for the eigenvalue solver

    unsigned int coeff_for_ncv; //!< eps_ncv = number_of_requested_ev * coeff_for_ncv

    double eigen_solver_tolerance; //!< tolerance for eigenvalue solver [Ha]

    std::string output_type; //!< output type 

    unsigned int max_iteration_number; //!< maximum number of iterations for the eigenvalue solver

    double spectrum_shift; //!< shift of spectrum ised in matrix assembly[eV]

    bool  consider_strain; //!< apply strain effect to the EFA Hamiltonian;

    bool  consider_potential; //!< apply strain effect to the EFA Hamiltonian;

    double disturb_arnoldi; //!< small disturb of the Hamiltonian matrix;

    bool Dirichlet_bc_everywhere;//!< apply dirichlet boundary conditions at all boundaries


    bool solve_ev_problem_twice;//!< if true, calculate the first eigenvalue only and then run again
 
  };



  struct eigen_energy
  {
    double energy; //!< eigen energy [eV]
    unsigned int global_number; //< eigen vector
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
    \param equation_systems reference to the "global" equation systems object
    \param opt  parameters  of the model
    \param problem_nam name of the problem that will be asigned to a new system
   
  */
  EnvelopFunctionApprox(EquationSystems&  equation_systems, std::string& problem_name, options& opt);


  //!destructor
  ~EnvelopFunctionApprox();


  //!computes Hamiltonian and S matrix
  void calculate_Hamiltonian_and_S(void); 



  //!sets material data
  void set_material_parameters(std::map<unsigned int, EFAbulkHamiltonian*>&  bulkHamiltonian); 


  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
   
  */
  void solve_eigen_value_problem(unsigned int ev_number );

  //!writes on disk the eigenfunction
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_eigen_function(unsigned int state_number,  std::string& filename);


  //!writes on disk the probability function \f$  \sum_i |\psi({\bf r})|^2   \f$
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_probability_function(unsigned int state_number,  std::string& filename);

  //! assigned mesh_data_objects
  /*!
    \param  mesh_data_in mesh data that contains material information and Dirichlet conditions
   */
  void assign_mesh_data(MeshData& mesh_data_in);

  


  //! define Dirichlet boundary condition
  /*!
    \param dirichlet_nodes_input nodes where \f$ \psi({\bf r}) = 0 \f$
  */
  void define_diriclet_nodes(std::vector<unsigned int>&  dirichlet_nodes_input);


  //! Passes a pointer of a Macrostrain object
  /*!
    \param strain_in a pointer to Macrostrain object
  */
  void define_strain_data( Macrostrain*  strain_in);



  //! Passes a pointer of a Drift-Diffusion object
  void define_Poisson_data( DriftDiffusion*  drift_in);


  //! returns  a pointer to the EquationSystems object
  EquationSystems* get_equation_systems();

  //! returns  a reference to the options
  const EnvelopFunctionApprox::options& get_options()  const;

  //! returns  a reference to the material_numbers object
  const std::vector<unsigned int>& get_material_numbers() const;

  //! returns a reference to solutions
  const std::vector<eigen_propblem_solution> get_solution() const;
 
  //! set spectrum shift  
  void set_spectrum_shift(double energy);

  //! returns conduction band minima for holes and valence band maximum for holes 
  double get_band_edge() const;

 private:



  //!Apply Dirichlet boundary conditions to all boundary points!
  void apply_diriclet_bc_at_all_boundaries();
  

  //!pointer to mesh of the equation systems
  Mesh* mesh;

  options opt;

  //!pointer to meshdata of the equation systems
  MeshData*  meshdata;

  EquationSystems* es;

  string system_name;

  //pointer to a drift-diffusion object that is used to get potential data 
  DriftDiffusion* poisson_equation;

  //pointer to the macrostrain object that is used to get strain data 
  Macrostrain* strain;

  //!system that we add to the equation systems
  LinearImplicitSystem* system;

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

  //!map that contains pointers to bulk Hamiltoninas
  std::map<unsigned int, EFAbulkHamiltonian*>  bulkHamiltonian;

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
  /*!
    \param number_of_ev number of eigen functions to read
  */
  void read_SLEPC_solution(unsigned int number_of_ev);

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
  void make_new_dofs(void);

  //!number of independent dofs
  int number_of_new_dofs;


  //!total number of dofs 
  int number_of_all_dofs;


  //! Hartree energy in eV
  static const double Hartree = 27.2113961;


  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons(eigen_propblem_solution state1, eigen_propblem_solution state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes(eigen_propblem_solution state1, eigen_propblem_solution state2);

  
  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons1(eigen_energy state1, eigen_energy state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes1(eigen_energy state1, eigen_energy state2);
  

  //! solutions of the eigenvalue problem
  std::vector<eigen_propblem_solution> solution;
  
};
#endif
