// $Id$

#ifndef _TIBERPETSCNONLINEARSOLVER_H_
#define _TIBERPETSCNONLINEARSOLVER_H_

#include "PetscRuntimeError.h"


// Libmesh includes
#include "nonlinear_solver.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"


#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscsnes.h>
}
#else
# include <petscsnes.h>
#endif



/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p NonlinearSolver<>
 *
 */
template <typename T>
class TiberPetscNonlinearSolver : public NonlinearSolver<T>
{
  public:

    //! Constructor. Initializes Petsc data structures
    TiberPetscNonlinearSolver(void) throw (PetscRuntimeError);

     //! Destructor.
    ~TiberPetscNonlinearSolver(void) throw (PetscRuntimeError);

    //! Release all memory and clear data structures.
    virtual void clear(void) throw (PetscRuntimeError);

    //! Initialize data structures if not done so already.
    virtual void init(void) throw (PetscRuntimeError);

    //!Call the Petsc solver.
    virtual std::pair<unsigned int, Real> solve(
        SparseMatrix<T>& jacobian,
        NumericVector<T>& solution,
        NumericVector<T>& residual,
        const double rtol,
        const unsigned int iter)
      throw (PetscRuntimeError, KSPDivergedError, SNESDivergedError);

    //! Set the options for the nonlinear solver
    void set_snes_options(double rtol = 1e-9, unsigned int max_it = 10);

    //! Set the options for the nonlinear solver
    void set_snes_options(double rtol, double atol, double stol,
        unsigned int max_it = 10);

    //! Set the options for the nonlinear solver
    /*!
     * \param ls_type the line search type:
     *  \li 1 = no linesearch (standard newton)
     *  \li 2 = quadratic line search
     *  \li 3 = cubic line search
     * \param ls_maxstep the maximum L2 norm for the line search step. This
     * is only used for \c ls_type 2 and 3
     */
    void set_snes_ls_options(int ls_type = 3, double ls_maxstep = 1e3);

    //! Set the options for the linear solver
    void set_ksp_options(double rtol = 1e-6, unsigned int max_it = 1000);

    //! Set the options for the linear solver
    void set_ksp_options(double rtol, double atol, unsigned int max_it = 1000);

    //! Set the linear solver (KSP) type
    void set_ksp_type(KSPType ksp_type);

    //! Set the preconditioner type
    /*!
     * Multiplicative composite preconditioning is supported by setting
     * \c PCCOMPOSITE as preconditioner type. It will use \c PCJACOBI
     * and \c PCILU
     */
    void set_pc_type(PCType pc_type);



  private:

    double _nonlinear_rtol;
    double _nonlinear_atol;
    double _nonlinear_stol;
    int _nonlinear_max_it;

    //! Is used to intercept spurious solver failures in the first iteration
    double _emergency_fnorm;
    
    double _linear_rtol;
    double _linear_atol;
    int _linear_max_it;

    int _ls_type;
    double _ls_maxstep;
    KSPType _ksp_type;
    PCType _pc_type;

    /**
     * Nonlinear solver context
     */
    SNES _snes;

    static void _checkerr(int errorcode) throw (PetscRuntimeError);

};


//
// inline methods
// 

template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::_checkerr(int errorcode)
  throw (PetscRuntimeError)
{
  if (errorcode != 0)
    throw(PetscRuntimeError(errorcode));
}

template <typename T>
inline
TiberPetscNonlinearSolver<T>::~TiberPetscNonlinearSolver(void)
  throw (PetscRuntimeError) 
{
  this->clear();
}

template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::set_snes_options(double rtol, unsigned int max_it)
{
  set_snes_options(rtol, 1e-15, 1e-6, max_it);
}


template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::set_ksp_options(double rtol,
    unsigned int max_it)
{
  set_ksp_options(rtol, 1e-50, max_it);
}

template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::set_ksp_options(double rtol, double atol,
    unsigned int max_it)
{
  _linear_rtol = rtol;
  _linear_atol = atol;
  _linear_max_it = max_it;
}


template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::set_ksp_type(KSPType ksp_type)
{
  _ksp_type = ksp_type;
}

template <typename T>
inline
void
TiberPetscNonlinearSolver<T>::set_pc_type(PCType pc_type)
{
  _pc_type = pc_type;
}



#endif // TIBERPETSCNONLINEARSOLVER_H_
