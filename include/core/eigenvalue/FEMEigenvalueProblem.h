// $Id$

#ifndef _FEMEIGENVALUEPROBLEM_H_
#define _FEMEIGENVALUEPROBLEM_H_

#include "EigenvalueProblem.h"



// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

#include <vector>
#include <complex>
#include <set>

class EquationSystems;
class MeshBase;
class LinearImplicitSystem;

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
   
    unsigned int number_of_eigenstates; //!< number of eigenstates to be calculated
    
    Method discretization_method;       //!< box integration or finite element
   
    bool Dirichlet_bc_everywhere;       //!< apply dirichlet boundary conditions at all boundaries

    std::string preconditioner;         //!< preconditioner name

    bool periodicity[3];                //!< periodic boundary conditions
  
    std::string spectral_trans;         //!< spectral transformation

    double spectrum_shift;              //!< Spectrum shift 

    std::string st_ksp_type;            //!< Liner system solution method 
    
    std::string strategy;               //<! matlab (algorithm used in Matlab) or general (recommended by SLEPC)

    bool monitor;   //<! activates convergence monitor if true

    double spectrum_inversion_tolerance; //<! tolerance for spectrum inversion

    bool dump_on_file;
    
  };


  FEMEigenvalueProblem(const ModelOptions& options);

  virtual ~FEMEigenvalueProblem() {};

  //! get H and S
  virtual int get_H_dim() const; 
  
  virtual int get_H_nnz() const;

  virtual void get_H_csr(std::vector<Complex>& A,std::vector<int>& JA,std::vector<int>& IA) const;

  virtual void get_S_csr(std::vector<Complex>& A, std::vector<int>& JA,std::vector<int>& IA) const;

  virtual void print_H(const std::string& outpath) const;
 

 protected:

  virtual void do_init();

  virtual void do_solve() {};

  virtual void parse_options();

  virtual void do_copy_H_to_solver();
  
  virtual void do_copy_S_to_solver();

  //!options of any eigensolver problem
  options solver_opt;


  //!dimension of the system
  short dim;


  //!number of independent dofs
  int number_of_new_dofs;


  //!total number of dofs 
  int number_of_all_dofs;


  //!vector: each element contains information about dof
  std::vector<FEMEigenvalueProblem::dof_new> new_dofs;

  //! my copy of DofMap::_dof_constraints (I have to recalculate it because it is private)
  DofConstraints 	my_dof_constraints;



  //!pointer to mesh of the equation systems
  MeshBase* mesh;


   //!pointer the equation systems object
  EquationSystems* es;
  
  //!system that we add to the equation systems
  LinearImplicitSystem* system;


  //!diriclet DOFS
  std::set<unsigned int>  dirichlet_dofs;

  //!list of periodic nodes
  //dim node list's: each contains list of nodes that periodic b.c
  //must be applied to
  std::vector<std::vector <const Node*>>  nodes_periodic; 


  //!creates new_dofs vector
  void make_new_dofs(void);

  //!updates my_dof_constraints 
  void make_constraints(void);

  //! checks if element lies on boundary
  bool element_on_boundary(const Elem* element);
    
  //!creates dirichlet dofs
  void create_dirichlet_dofs(void);
  
  //!Apply Dirichlet boundary conditions to all boundary points!
  void apply_dirichlet_at_all_boundaries();
  
  //! Apply periodic boundary conditions
  void apply_periodic_bc();

  //! create list of nodes that lies at the periodic boundary
  void make_nodes_periodic();

  //! Apply bc to eigenvalue problem
  void apply_bc();
                                                                 
  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spectrum shift 
    \param haveS decides whether it is a generalize eigenproblem
  */
  void solve_eigen_value_problem(unsigned int ev_number, double spectrum_shift = 0.0);


  //!simulation domain boundary
  double min_coord[3];
  
  //!simulation domain boundary
  double max_coord[3];


  //! Get periodicity vector of a constrained DOF
  Point get_periodicity_vector(unsigned int dof) const;


 private:

  //!size of the Hamiltonian matrix
  unsigned int _hamiltonian_size;

  std::map<unsigned int, Point> _constrained_dof_periodicity;


};


inline
Point
FEMEigenvalueProblem::get_periodicity_vector(unsigned int dof) const
{
  Point p(0);
  std::map<unsigned int, Point>::const_iterator it(
      _constrained_dof_periodicity.find(dof));
  if (it != _constrained_dof_periodicity.end())
    p = it->second;

  return p;
}


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
