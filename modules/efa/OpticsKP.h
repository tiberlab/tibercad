#ifndef _OpticsKP_H_
#define _OpticsKP_H_


#include "tibercad/physics/optics/Optics.h"

#include "EnvelopFunctionApprox.h"
#include "tibercad/kintegration/KspaceIntegration.h"
#include "tibercad/kintegration/KspaceIntegrationTemplate.h"




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
  libMesh::LinearImplicitSystem* system{nullptr};

  std::vector<unsigned int> psivar;
  

  //!pointer to the EFA for initial states to access its class members
  EnvelopFunctionApprox* initial_state_model{nullptr};


  //!pointer to the EFA for initial states to access its class members 
  EnvelopFunctionApprox* final_state_model{nullptr};


  //!map that contains pointers to bulk Hamiltoninas
  std::map<unsigned int, KPbulkHamiltonian*>  bulkHamiltonian;

  
  //!pointer to the real part of  Px matrix
  libMesh::SparseMatrix<Number>* Px_matr_real{nullptr};

  //!pointer to the imaginary part of  Px matrix
  libMesh::SparseMatrix<Number>* Px_matr_imag{nullptr};

  //!pointer to the real part of  Py matrix
  libMesh::SparseMatrix<Number>* Py_matr_real{nullptr};
 
  //!pointer to the imaginary part of  Py matrix
  libMesh::SparseMatrix<Number>* Py_matr_imag{nullptr};
 
  //!pointer to the real part of  Pz matrix
  libMesh::SparseMatrix<Number>* Pz_matr_real{nullptr};

  //!pointer to the imaginary part of  Pz matrix
  libMesh::SparseMatrix<Number>* Pz_matr_imag{nullptr};
  
  /*!
   * \brief Whether or not to compute the envelope
   *        contribution to the optical matrix element
   *
   * As default, we include both the overlap and the momentum
   * matrix element of the envelope functions, as it is no
   * large computational overhead.
   */
  bool _compute_intraband_terms{true};




  //!calculate P-vector matrix element between states i and j
  /*!
    \param i initial state number
    \param j final state number
  */
  std::vector<libMesh::Complex>  calculate_matrix_element(unsigned int i, unsigned int j);






};

inline OpticsKP* OpticsKP::create(const ModelOptions& options)
{
  return (new OpticsKP(options));
}



 



#endif
