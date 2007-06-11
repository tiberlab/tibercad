// $Id$


#ifndef _TIBERPETSCLINEARSOLVER_H_
#define _TIBERPETSCLINEARSOLVER_H_

#include "PetscRuntimeError.h"

// Libmesh includes
#include "linear_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"



#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petscksp.h>
}
#else
# include <petscversion.h>
# include <petscksp.h>
#endif



//! The TiberCAD linear solver interface to PETSc
/*!
 * This class provides the TiberCAD interface to PETSc
 * iterative solvers.
 * It is derived from the libmesh LinearSolver class
 *
 */
template <typename T>
class TiberPetscLinearSolver : public LinearSolver<T>
{

  public:

    //!  Constructor. Initializes Petsc data structures
    TiberPetscLinearSolver(void);


    //! Destructor.
    virtual ~TiberPetscLinearSolver(void);


    //! Release all memory and clear data structures.
    virtual void clear(void);


    //! Initialize data structures if not done so already.
    virtual void init(void);


    //! Set the options for the linear solver
    void set_ksp_options(double rtol = 1e-6, unsigned int max_it = 1000);


    //! Set the options for the linear solver
    void set_ksp_options(double rtol, double atol, unsigned int max_it = 1000);


    //! Call the Petsc solver.
    /*!
     * It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */    
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<T>  &matrix_in,
          NumericVector<T> &solution_in,
          NumericVector<T> &rhs_in,
          const double tol,
          const unsigned int m_its)
      {
        return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
      }


    //! Call the linear solver specifying explicitly the preconditioner matrix
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<T>  &matrix,
          SparseMatrix<T>  &preconditioner,
          NumericVector<T> &solution,
          NumericVector<T> &rhs,
          const double tol,
          const unsigned int m_its);


    /*!
     * \brief Fills the input vector with the sequence of residual norms
     * from the latest iterative solve.
     */
    void get_residual_history(std::vector<double>& hist);


    //! Returns just the initial residual for the solve just completed
    double get_initial_residual();


  private:


    double _linear_rtol;
    double _linear_atol;
    int _linear_max_it;


    //! Preconditioner context
    PC _pc; 


    //! Krylov subspace context
    KSP _ksp;


    //! Set the KSP type
    void set_ksp_type(void);


    //! Set the KSP type
    void set_pc_type(void);


    //! Check PETSc error code
    static void _checkerr(int errorcode) throw (PetscRuntimeError);

};



//
// inline members
//


template <typename T>
inline
void
TiberPetscLinearSolver<T>::_checkerr(int errorcode) throw (PetscRuntimeError)
{
  if (errorcode != 0)
    throw(PetscRuntimeError(errorcode));
}





template <typename T>
inline
TiberPetscLinearSolver<T>::TiberPetscLinearSolver(void)
  : _linear_rtol(1e-6),
    _linear_atol(1e-50),
    _linear_max_it(500)
{
  if (libMesh::n_processors() == 1)
    this->_preconditioner_type = ILU_PRECOND;
  else
    this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
}



template <typename T>
inline
TiberPetscLinearSolver<T>::~TiberPetscLinearSolver(void)
{
  this->clear();
}



template <typename T>
inline
void
TiberPetscLinearSolver<T>::set_ksp_options(double rtol,
    unsigned int max_it)
{
  set_ksp_options(rtol, 1e-50, max_it);
}



template <typename T>
inline
void
TiberPetscLinearSolver<T>::set_ksp_options(double rtol, double atol,
    unsigned int max_it)
{
  _linear_rtol = rtol;
  _linear_atol = atol;
  _linear_max_it = max_it;
}




#endif // _TIBERPETSCLINEARSOLVER_H_
