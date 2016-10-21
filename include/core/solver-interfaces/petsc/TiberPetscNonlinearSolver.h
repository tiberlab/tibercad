// $Id$

#ifndef _TIBERPETSCNONLINEARSOLVER_H_
#define _TIBERPETSCNONLINEARSOLVER_H_

#include "TiberNonlinearSolver.h"
#include "PetscRuntimeError.h"


// Libmesh includes
#include "petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
# include <petscsnes.h>
EXTERN_C_FOR_PETSC_END


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
    TiberPetscNonlinearSolver(sys_type& s);

     //! Destructor.
    ~TiberPetscNonlinearSolver(void);

    //! Release all memory and clear data structures.
    virtual void clear(void);

    //! Initialize data structures if not done so already.
    virtual void init(const char* name = NULL);
    //virtual void init(void);
    
    //!Call the Petsc solver.
    virtual std::pair<unsigned int, Real> solve(
        SparseMatrix<double>& jacobian,
        libMesh::NumericVector<double>& solution,
        libMesh::NumericVector<double>& residual);

    //! Get the divergence tolerance
    double get_divergence_tol(void) const;

    //! Set or get the old gnorm
    double& old_gnorm(void);


  protected:

    /*! \copydoc TiberNonlinearSolver::parse_options() */
    virtual void parse_options(const ModelOptions& options);
    
    virtual int get_total_linear_iterations ();
    virtual unsigned get_current_nonlinear_iteration_number () const;


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
    std::string _ksp_type;

    //! The PC type
    std::string _pc_type;

    //! The linear relative tolerance
    double _linear_rtol;

    //! The linear absolute tolerance
    double _linear_atol;

    //! The maximum number of linear iterations
    int _linear_max_it;


    //! Nonlinear solver context
    SNES _snes;

};


//
// inline methods
// 


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


inline
int
TiberPetscNonlinearSolver::get_total_linear_iterations(void)
{
  return 0;
}

inline
unsigned
TiberPetscNonlinearSolver::get_current_nonlinear_iteration_number(void) const
{
  return 0;
}


#endif // TIBERPETSCNONLINEARSOLVER_H_
