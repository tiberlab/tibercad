#ifndef _OpticsKP_H_
#define _OpticsKP_H_

#include  "EnvelopFunctionApprox.h"
#include "KPbulkHamiltonian.h"




//!Class that calculates optical matrix elements in k.p formalism
class OpticsKP
{
 public:


  typedef std::complex<double> Complex;

  //! constructor
  OpticsKP();


  //! constructor
  /*!
    \param initial_state_model intial states of the optical transition
    \param final_states_model intial states of the optical transition
   */
  OpticsKP(const EnvelopFunctionApprox* initial_state_model, const EnvelopFunctionApprox* final_state_model, EquationSystems* es, std::string& system_name);


  ~OpticsKP();
  
  void set_initial_eigen_states(const std::vector < unsigned int >& initial_eigen_state_numbers	); 

  void set_final_eigen_states(const std::vector < unsigned int >& final_eigen_state_numbers	); 

  //! calculate momentum matrix elements 
  /*!
    \param P_matrix matrix that contains momentum matrix elements
    P_matrix[i][j][k]: i - component (x,y,z); j - initial state; k - final state
  */
  void get_P_matrix_elements (std::vector< std::vector <std::vector <  Complex  >  >  > &  P_matrix);


  //! set kp material model for optics calculation
  /*!
    Has to 8x8 model. Parameters may be different from those used in Schroedinger equation.
   */ 
  void set_material_parameters(std::map<unsigned int, KPbulkHamiltonian*>&  bulkHamiltonian); 

  //! calculate Px, Py and Pz matrixes 
  void calculate_matrix(void);


  //!calculate averaged value of the electrochemical potential <\psi|\mu|psi>
  /*!
    \param  i number of state
    \param  kind  1: the initial state; 2: the final state
   */
  double calculate_fermi_averaged(unsigned int i, short kind);


 private:

  //!equation systems object
  EquationSystems* es ;
  

  std::vector<unsigned int> _initial_eigen_state_numbers;


  std::vector<unsigned int> _final_eigen_state_numbers;


  const EnvelopFunctionApprox* initial_state_model;



  const EnvelopFunctionApprox* final_state_model;

 

  //!map that contains pointers to bulk Hamiltoninas
  std::map<unsigned int, KPbulkHamiltonian*>  bulkHamiltonian;

  
  //!pointer to the real part of  Px matrix
  SparseMatrix<Number>* Px_matr_real;

  //!pointer to the imaginary part of  Px matrix
  SparseMatrix<Number>* Px_matr_imag;

  //!pointer to the real part of  Py matrix
  SparseMatrix<Number>* Py_matr_real;
 
  //!pointer to the imaginary part of  Py matrix
  SparseMatrix<Number>* Py_matr_imag;
 
  //!pointer to the real part of  Pz matrix
  SparseMatrix<Number>* Pz_matr_real;

  //!pointer to the imaginary part of  Pz matrix
  SparseMatrix<Number>* Pz_matr_imag;

  
  //!name of the system
  std::string system_name;


  //!system that we add to the equation systems
  LinearImplicitSystem* system;


  //!my Jacobian because I calculate everything in atomic units
  double my_Jacobian; 

  //!calculate P-vector matrix element between states i and j
  /*!
    \param i initial state number
    \param j final state number
  */
  std::vector<Complex> calculate_matrix_element(unsigned int i, unsigned int j);

 


};

#endif
