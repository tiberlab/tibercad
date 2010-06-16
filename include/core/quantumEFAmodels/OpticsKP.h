#ifndef _OpticsKP_H_
#define _OpticsKP_H_


#include "Optics.h"

#include "EnvelopFunctionApprox.h"
class  KPbulkHamiltonian;




//!Class that calculates optical matrix elements in k.p formalism
class OpticsKP: public Optics
{
 public:

  //!control options
  enum JobKind
  {
    MATREL = 0, //!< calculate matrix elements of momentum matrix
    BULKMATREL = 1, //!<  calculate matrix elements of momentum matrix for bulk simulation
  };




  typedef std::complex<double> Complex;

  //! constructor
  OpticsKP(const ModelOptions& options);

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
    create_physical_model(const ModelOptions& options,
        const Material* mat) const throw (ModelErrorException);
    
  
  //!Here returns NULL
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);



  // void set_material_parameters(std::map<unsigned int, KPbulkHamiltonian*>&  bulkHamiltonian); 

  //! calculate Px, Py and Pz matrixes 
  void calculate_matrix(void);



  //! calculate Px, Py and Pz matrixes for bulk 
  void calculate_matrix_bulk(void);



 

 

  static OpticsKP* create(const ModelOptions& options);
  
  //!calculate spectrum 
  /*!
    \f$
    
    P(\hbar \omega) = \sum_{i,j} \frac{1}{2\pi^2}  \frac{\omega^2_{ij} e^2 }{m^2 c^3}  |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
    \frac{\Gamma/2} {(\hbar \omega_{ij} - \hbar \omega)^2 + (\Gamma/2)^2} d\Omega
    \f$
      
    \param Energy energy grid [eV]
    \param spectrum calculated spectrum (atomic units)
    \param Gamma broadering parameter [eV]
    \param polariz polarization vector of a linearly polarized light (must be a normalized one, \f$ |{\bf e}| = 1 \f$)
    


  */
  virtual void calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
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

  std::vector<unsigned int> psivar;
  

  //! Mesh for spectrum [eV];
  Mesh* _energy_mesh;

  //!defines which job has to be done
  JobKind job;


 protected:

  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


 //! plot the  spectrum results for the single k-point case
  virtual void do_plot();


};

inline OpticsKP* OpticsKP::create(const ModelOptions& options)
{
  return (new OpticsKP(options));
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
