#ifndef _OpticsKP_H_
#define _OpticsKP_H_


#include "Optics.h"

#include "EnvelopFunctionApprox.h"
#include "KspaceIntegration.h"
#include "KspaceIntegrationTemplate.h"

class  KPbulkHamiltonian;
class  Device;

 


//!Class that calculates optical matrix elements in k.p formalism
class OpticsKP: public Optics
{
 public:


  //! constructor
  OpticsKP(const ModelOptions& options);

  ~OpticsKP();

  static OpticsKP* create(const ModelOptions& options);
 
 


  virtual PhysicalModel*
    create_bulk_model(const ModelOptions& options,
        const Material* mat) const;
    
  

 protected:

  virtual void 	do_init(void);

  //! calculate Px, Py and Pz matrixes 
  virtual void do_assemble(const ModelOptions& opts);

  //! calculate Px, Py and Pz matrixes for bulk 
  virtual void calculate_matrix_bulk(void);

  virtual void do_compute_matrix_elements(void);

 private:

  //!system that we add to the equation systems
  LinearImplicitSystem* system;

  std::vector<unsigned int> psivar;
  

  //!pointer to the EFA for initial states to access its class members
  EnvelopFunctionApprox* initial_state_model;


  //!pointer to the EFA for initial states to access its class members 
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
  

  //!calculate P-vector matrix element between states i and j
  /*!
    \param i initial state number
    \param j final state number
  */
  std::vector<Complex>  calculate_matrix_element(unsigned int i, unsigned int j);



  //! creates kp material model for optics calculation
  /*!
    Has to be 8x8 model. Parameters may be different from those used in Schroedinger equation.
  */ 




};

inline OpticsKP* OpticsKP::create(const ModelOptions& options)
{
  return (new OpticsKP(options));
}



 



#endif
