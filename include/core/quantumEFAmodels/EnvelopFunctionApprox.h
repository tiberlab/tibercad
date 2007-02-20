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
#include "Macrostrain.h"
#include "DriftDiffusion.h"
#include "SimulationInterface.h"

class EnvelopFunctionApprox  : public SimulationInterface
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


    bool periodicity[3];    //!< periodic boundary conditions

    std::string solver;     //!< solver type

    std::string mpi_command_line; //!< something like:  mpirun -np 5 -machinefile machines

    std::string solver_command_line; //!< additional line to pass parameters for the eigenvalue solver

    unsigned int coeff_for_ncv; //!< eps_ncv = number_of_requested_ev * coeff_for_ncv

    double eigen_solver_tolerance; //!< tolerance for eigenvalue solver [Ha]

    std::string output_type; //!< output type 

    unsigned int max_iteration_number; //!< maximum number of iterations for the eigenvalue solver

    double spectrum_shift;    //!< shift of spectrum ised in matrix assembly[eV]

    bool  consider_strain;    //!< apply strain effect to the EFA Hamiltonian;

    bool  consider_potential; //!< apply strain effect to the EFA Hamiltonian;

    bool estimate_spectrum_shift; //!< calculate spectrum shift from band edges;
  

    bool Dirichlet_bc_everywhere;//!< apply dirichlet boundary conditions at all boundaries


    bool solve_ev_problem_twice;//!< if true, calculate the first eigenvalue only and then run again


    bool convergent_density;//!< if true, the number of eigenstates will be increased to reach the tolerance

 
    unsigned int initial_eigenstates_number; //!< initial number of eigenstates that is used in an iterative calculation of the density


    double relative_density_tolerance; //!< stops itarations if \f$ \rho_i / \rho_{i+1} < \varepsilon    \f$, where \f$ \rho \f$ is the                                              total density 
 
    double eigen_number_increase_factor; //!< to increase number of eigenstates for the next iteration 


    bool log_output; //!< to do a lot of output on screen


    unsigned int number_of_eigenstates; //!<number of eigenstates to be calculated

    

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
    double Fermi_energy; //< electro-chemical potential [eV] \f$ \langle \psi |\mu({\bf r} | \psi \rangle  \f$
  };


  //! data structure that contain information about if the dof is independent or not
  struct dof_new
  {
    bool independent;  //!< true if it is and independent dof
    unsigned int new_number;  //!< new number in the independent dofs list
  };


 

  //!constructor
  EnvelopFunctionApprox(void);


  //!destructor
  ~EnvelopFunctionApprox();


  //!computes Hamiltonian and S matrix
  void calculate_Hamiltonian_and_S(void); 



  

  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spetrum shift [eV]
  */
  void solve_eigen_value_problem(unsigned int ev_number, double spectrum_shift = 0.0 );

  //!writes on disk the eigenfunction
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_eigen_function(unsigned int state_number,  const std::string& filename);


  //!writes on disk the probability function \f$  \sum_i |\psi({\bf r})|^2   \f$
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_probability_function(unsigned int state_number,  const std::string& filename);


  

  //! returns a reference to solutions
  const std::vector<eigen_propblem_solution>& get_solution() const;
 
 
  //! returns conduction band minima for holes and valence band maximum for holes 
  double get_band_edge() const;


  
  //!calculate averaged value of the electrochemical potential <\psi|\mu|psi>
  /*!
    \param  i number of state
   
   */
  double calculate_fermi_averaged(unsigned int i);




  
  //!apply k Block vector to all the Hamiltonians
  /*!
    \param k k-vector in atomic units
  */
  void apply_k_vector(const double k[3]);



  //! claculate total density
  /*!
    \f$ \rho = \sum_i F_{\rm{Fermi}}(E_i) \f$
    \param T temperature [K];
  */
  double get_integrated_probability(double T);



  //!obtain convergent density
  /*!
    \param T temperature [K]
    \param cell_data  if true, then cell data is calculated (default); if false, the nodal data ic calculated 
  */
  std::vector<double>  calculate_convergent_density(double T, bool cell_data = true);
  


   //! calculate nodal or cell density  (in atomic units) for a single \f$ {\bf k}_{\|}\f$ vector.  
  /*!
    The nodal density reads: \f$ \rho({\bf r}) = \sum_i   |\psi_i({\bf r})|^2 F_{fermi}(E_i) \f$ 
    The cell density reads:  \f$ \rho = \frac{1}{\Omega_0} \int_{\Omega_0}  \sum_i   |\psi_i({\bf r})|^2 \, dV  F_{fermi}(E_i), \f$
    where \f$ \Omega_0 \f$ is the element volume.

    \param T temperature [K]
    \param cell_data  if true, then cell data is calculated (default); if false, the nodal data ic calculated 

  */
  std::vector<double>  calculate_density(double T, bool cell_data = true);


  //! sets opt.initial_eigestates_number
  void set_initial_eigenstates_number(unsigned int n);



  //! returns number of active cells
  unsigned int get_number_of_active_cells();


  static  EnvelopFunctionApprox* create();


  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException);
    
   
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);

 private:

  //!pointer to the device object
  static  Device* _device;


  //!Apply Dirichlet boundary conditions to all boundary points!
  void apply_diriclet_bc_at_all_boundaries();
  

  //!pointer to mesh of the equation systems
  Mesh* mesh;

  options opt;

  //!pointer to meshdata of the equation systems
  MeshData*  meshdata;

  //!pointer the equation systems object
  EquationSystems* es;

  
  std::string system_name;

  //!pointer to a drift-diffusion object that is used to get potential data 
  DriftDiffusion* poisson_equation;

  //!pointer to the macrostrain object that is used to get strain data 
  Macrostrain* strain;

  //!system that we add to the equation systems
  LinearImplicitSystem* system;

  //!diriclet nodes vector
  // std::vector<unsigned int>  dirichlet_nodes;

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
  //std::vector<unsigned int> material_of_elem;


  //!create material_number vector
  //void assemble_material_list(void);


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
 
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want 
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize eigenfunctions
  5) calculate fermi energy for each state
 
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
  


  //! my copy of DofMap::_dof_constraints (I have to recalculate it because it is private)
  DofConstraints 	my_dof_constraints;

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
  
  //! cheks if element lies on boundary
  bool element_on_boundary(const Elem* element);
  


  //! calculates the norm of the eigenstate \f$ \sqrt {| \langle \psi|\psi \rangle |} \f$
  /*!
    \param state_number number of the eigenstate
  */
  double eigenstate_norm(unsigned int state_number);


 
  
  //! calculate density without \f$ | \psi_i (r) |^2 \f$
  /*!
    \param i number of the eigenstate
   */
  std::vector<double> calculate_prob_function(unsigned int i);


  //! calculate density without \f$ \frac{1}{\Sigma_0} \int_{\Sigma_0} | \psi_i (r) |^2 /,dV \f$
  /*!
    \param i number of the eigenstate
  */
  std::vector<double> calculate_cell_prob_function(unsigned int i);



  //!Calculates Fermi Dirac probability
  /*!
    For electrons:  \f$ p = \frac{1}{1 + \exp (\frac{E - \mu}{kT})}     \f$
    For holes:      \f$ p = 1 - \frac{1}{1 + \exp (\frac{E - \mu}{kT})} \f$
    \param Energy   state energy [eV]
    \param Fermi_energy  Fermi energy [eV]
    \param temperature   temperature [K] 
   */
  double Fermi_statistics_probability(double Energy, double Fermi_energy, double temperature);


  //!Calculate number of bands in the Hamiltonian 
  short calculate_number_of_bands(void) const;



 

 protected:


  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


 



};
//-------------------------------------------------------------------
inline double EnvelopFunctionApprox::Fermi_statistics_probability(double Energy, double Fermi_energy, double Temperature)
{
  

  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;

  if (exp_arg > 20) 
    el_fermi = 0.0;
  else
    el_fermi = 1.0/(  1 +  std::exp(exp_arg)  );


  

  if (opt.particle == "el")
    return(el_fermi);
  else
    return(1.0 - el_fermi);


}

//-------------------------------------------------------------------
inline bool EnvelopFunctionApprox::element_on_boundary(const Elem* element)
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
  
};
//---------------------------------------------------------

inline void EnvelopFunctionApprox::set_initial_eigenstates_number(unsigned int n)
{
  opt.initial_eigenstates_number = n;
}


//---------------------------------------------------------

inline EnvelopFunctionApprox*  EnvelopFunctionApprox::create()
{
  return (new EnvelopFunctionApprox );
}
#endif
