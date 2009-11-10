// $Id$

#ifndef _TIBERPETSCNONLINEARSOLVER_H_
#define _TIBERPETSCNONLINEARSOLVER_H_

#include "PetscRuntimeError.h"


// Libmesh includes
#include "TiberNonlinearSolver.h"


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
class TiberPetscNonlinearSolver : public TiberNonlinearSolver
{
  public:

    //! Constructor. Initializes Petsc data structures
    TiberPetscNonlinearSolver(void);

     //! Destructor.
    ~TiberPetscNonlinearSolver(void);

    //! Release all memory and clear data structures.
    virtual void clear(void);

    //! Initialize data structures if not done so already.
    virtual void init(void);

    //!Call the Petsc solver.
    virtual std::pair<unsigned int, Real> solve(
        SparseMatrix<double>& jacobian,
        NumericVector<double>& solution,
        NumericVector<double>& residual);

    //! Get the divergence tolerance
    double get_divergence_tol(void) const;

    //! Set or get the old gnorm
    double& old_gnorm(void);


  protected:

    /*! \copydoc TiberNonlinearSolver::parse_options() */
    virtual void parse_options(const ModelOptions& options);
    

  private:

    //! Is used to intercept spurious solver failures in the first iteration
    double _emergency_fnorm;
    
    //! The line search type
    int _ls_type;

    //! The maximum linesearch step
    double _ls_maxstep;

    //! The previous gnorm for divergence test
    double _old_gnorm;

    //! The divergence tolerance
    /*!
     * If gnorm > _divergence_tol * _old_gnorm we assume divergence of
     * algorithm
     */
    double _divergence_tol;


    //! The KSP type
    KSPType _ksp_type;

    //! The PC type
    PCType _pc_type;

    //! Nonlinear solver context
    SNES _snes;

    static void _checkerr(int errorcode);

};


//
// inline methods
// 

inline
void
TiberPetscNonlinearSolver::_checkerr(int errorcode)
{
  if (errorcode != 0)
    throw(PetscRuntimeError(errorcode));
}

inline
TiberPetscNonlinearSolver::~TiberPetscNonlinearSolver(void)
{
  this->clear();
}



inline
double
TiberPetscNonlinearSolver::get_divergence_tol(void) const
{
  return _divergence_tol;
}



inline
double&
TiberPetscNonlinearSolver::old_gnorm(void)
{
  return _old_gnorm;
}




#endif // TIBERPETSCNONLINEARSOLVER_H_
