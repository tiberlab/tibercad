#ifndef _MAXWELLEQUATIONS_H_
#define _MAXWELLEQUATIONS_H_
#include "Device.h"
#include "EigenvalueProblem.h"




//!Class to solve Maxwell equations
class MaxwellEquations : public EigenvalueProblem
{
 public:

  struct eigen_value
  {
    double k_squared; //!< \f$ \frac{\omega^2}{c^2} \f$
    unsigned int global_number; //< eigen vector number
  };



  //! data structure that contains eigenvalue and eigenvector
  struct eigen_problem_solution
  {
    double k_squared; // !< \f$ \frac{\omega^2}{c^2} \f$
    std::vector< std::complex<double>  > eigen_vector; //< eigen vector
  };
  


  //!constructor
  MaxwellEquations(void);
  
  //!destructor
  virtual ~MaxwellEquations(void) {};
  
  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException);
    
   
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);


  static  MaxwellEquations* create();




  struct options
  {
    double work_units; //!< work units for length [m] 
    double spectrum_shift ;
    bool scalar_approximation;
  
  };


 protected:


  //!in this class  outputs \f$ |{\bf E}|^2 \f$ for each eigenstate
  virtual void build_nodal_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);


  virtual void 	build_elemental_results (const std::set< std::string > &variables, 
					 std::vector< double > &results, std::vector< std::string > &legend){};


  //!in this class it outputs eigen photon energies
  virtual void 	build_integrated_quantities (const std::set< std::string > &names, std::vector< double > &values);

  //!in this class it builds description for eigen values
  virtual void 	build_integrated_quantities_description (const std::set< std::string > &names,
							 std::vector< std::string > &legend, 
							 std::vector< std::string > &description);
 


  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


 

 private:

  //!pointer to the device object
  static  Device* _device;
 
 

  //!name of the system
  std::string system_name;
  


  //!solver options
  options opt;


  //!computes Hamiltonian and S matrix
  void calculate_Hamiltonian_and_S(void);

  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spetrum shift [eV]
  */
  // void solve_eigen_value_problem(unsigned int ev_number, double spectrum_shift = 0.0 );
  

  //! Kronecker delta \f$ \delta_{ij} \f$ 
  inline int delta_Kronecker(int i, int j);
  
  
  //! Tensor product of two Levi-Civita sybmols:  \f$  \sum_{i = 1} ^3 e_{ijk}e_{imn} = \delta_{jm}\delta_{kn} - \delta{jn}\delta{km} \f$
  inline int LeviCivita_product(int j, int k , int m, int n);


  //!read SLEPc solutions
  /*!
 
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want 
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize
  
 
    \param number_of_ev number of eigen functions to read
  */
  virtual void read_SLEPC_solution(unsigned int number_of_ev);

  //! solutions of the eigenvalue problem
  std::vector<eigen_problem_solution> solution;

  //! spectrum shift is almost equal to the lowest eigenvalue
  virtual double get_new_spectrum_shift(void);


 
  void copy_S_matrix_to_solver(void);


  //! calculated \f$ {\bf E(r)}^2 \f$
  /*!
    \param mode_number number of eigenmode
    \param data  \f$ {\bf E(r)}^2
  */
  void prepare_field_mod_squared(const unsigned int mode_number, std::vector<double>& data);

  
  

  inline static bool compare_eigenstate(eigen_problem_solution state1, eigen_problem_solution state2);

 
  inline static bool compare_eigenvalue(eigen_value state1, eigen_value state2);
  

  //!calculates norm \f$ \left( \int  |{\bf A}|^2 dV  \right)^{1/2} \f$ 
  double  eigenstate_norm(unsigned int state_number);


  short number_of_field_components;


};



inline int MaxwellEquations::delta_Kronecker(int i, int j)
{
  return ((i==j) ? 1 : 0);
}


//------------------------------------------------------------------------//


inline int MaxwellEquations:: LeviCivita_product(int j, int k , int m, int n)
{
  return ( delta_Kronecker( j,m )*delta_Kronecker( k,n ) -  delta_Kronecker( j,n )*delta_Kronecker( k,m )  );
}


//------------------------------------------------------------------------//

inline MaxwellEquations* MaxwellEquations::create()
{
 return (new MaxwellEquations);
}

//------------------------------------------------------------------------//
inline  bool MaxwellEquations::compare_eigenstate(MaxwellEquations::eigen_problem_solution state1, 
						      MaxwellEquations::eigen_problem_solution state2)
{
  return(state1.k_squared < state2.k_squared );
}


//------------------------------------------------------------------------//
inline  bool MaxwellEquations::compare_eigenvalue( MaxwellEquations::eigen_value state1,  MaxwellEquations::eigen_value state2)
{
  return(state1.k_squared < state2.k_squared );
}

#endif
