#ifndef _OpticsKP_H_
#define _OpticsKP_H_


#include "SimulationInterface.h"

#include "EnvelopFunctionApprox.h"
class  KPbulkHamiltonian;




//!Class that calculates optical matrix elements in k.p formalism
class OpticsKP: public SimulationInterface
{
 public:


  typedef std::complex<double> Complex;

  //! constructor
  OpticsKP();

  ~OpticsKP();

 
 
  
 
  //! calculate momentum matrix elements 
  /*!
    \param P_matrix matrix that contains momentum matrix elements
    P_matrix[i][j][k]: i - component (x,y,z); j - initial state; k - final state
  */
  void calculate_P_matrix_elements (void);


  //! creates kp material model for optics calculation
  /*!
    Has to be 8x8 model. Parameters may be different from those used in Schroedinger equation.
  */ 

  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException);
    
  
  //!Here returns NULL
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);



  // void set_material_parameters(std::map<unsigned int, KPbulkHamiltonian*>&  bulkHamiltonian); 

  //! calculate Px, Py and Pz matrixes 
  void calculate_matrix(void);


 

  static OpticsKP* create();
  
  //!calculate spectrum 
  /*!
    \param Energy energy grid [eV]
    \param spectrum calculated spectrum (atomic units)
    \param Gamma broadering parameter [eV]
    \param polariz polarization vector of a linearly polarized light (must be a normalized one, \f$ |{\bf e}| = 1 \f$)
    \f$
    
    P(\hbar \omega) = \sum_{i,j} \omega |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
    \frac{\Gamma/2}{(E_i - E_j - \hbar \omega)^2 + (\Gamma/2)}
    \f$


  */
  void calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
                          std::map<const Elem*, double>& spectrum);

  //! get pointer to initial state simulation (to be  used by OpticalSpectrum)
  EnvelopFunctionApprox* get_initial_state_model();

  // get pointer to final state simulation (to be  used by OpticalSpectrum)
  EnvelopFunctionApprox* get_final_state_model();


 private:

  //!equation systems object
  EquationSystems* es ;
  

  //!pointer to the device object
  static  Device* _device;


  //!numbers of eigensates that are considered as intial states for optical transition
  std::vector<unsigned int> _initial_eigen_state_numbers;


  //!numbers of eigensates that are considered as final states for optical transition
  std::vector<unsigned int> _final_eigen_state_numbers;

  //!pointer to the eigenvalue solver for initial states
   EnvelopFunctionApprox* initial_state_model;


  //!pointer to the eigenvalue solver for final states
   EnvelopFunctionApprox* final_state_model;

 

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

  //!length scale from mesh units to atomic unit
  double length_scale;

  //!my Jacobian because I calculate everything in atomic units
  double my_Jacobian; 

  //!calculate P-vector matrix element between states i and j
  /*!
    \param i initial state number
    \param j final state number
  */
  std::vector<Complex> calculate_matrix_element(unsigned int i, unsigned int j);



  
  //!k-vector in atomic units
  double k_vector[3];

  //!  momentum matrix elements 
  /*!
    P_matrix[i][j][k]: i - component (x,y,z); j - initial state; k - final state
  */
  std::vector< std::vector <std::vector <  Complex  >  >  >   P_matrix;

 protected:

  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


};

inline OpticsKP* OpticsKP::create()
{
  return (new OpticsKP);
}

inline 
EnvelopFunctionApprox* OpticsKP::get_initial_state_model()
{
  return initial_state_model;

}

inline
EnvelopFunctionApprox* OpticsKP::get_final_state_model()
{
  return final_state_model;

}

 



#endif
