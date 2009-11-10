// $Id$


#ifndef _TIBERPETSCLINEARSOLVER_H_
#define _TIBERPETSCLINEARSOLVER_H_

#include "TiberLinearSolver.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"

// Libmesh includes
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
class TiberPetscLinearSolver : public TiberLinearSolver
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


    //! Solve the linear system
    virtual std::pair<unsigned int, Real>
      solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs);

    /*!
     * \brief Fills the input vector with the sequence of residual norms
     * from the latest iterative solve.
     */
    void get_residual_history(std::vector<double>& hist);


    //! Returns just the initial residual for the solve just completed
    double get_initial_residual();

    
    //! Returns KSP context
    KSP get_ksp(void);

    


  protected:

    /*!  \copydoc TiberLinearSolver::parse_options() */
    virtual void parse_options(const ModelOptions& options);


    //! Setup the textual and graphical convergence monitors
    void setup_monitors(void);



  private:


    //! Krylov subspace context
    KSP _ksp;

    //! The KSP type
    KSPType _ksp_type;

    //! The PC type
    PCType _pc_type;

    //! The graphical monitor
    PetscDrawLG _LG_monitor;

    //! Do we want monitor?
    bool _monitor;

    //! Do we want X monitor?
    bool _xmonitor;

    //! Tells if xmonitor is already opened
    bool _xmonitor_open;


    //! Check PETSc error code
    static void _checkerr(int errorcode);


    //! Check convergence
    std::pair<unsigned int, double> check_convergence(void);


};



//
// inline members
//


inline
void
TiberPetscLinearSolver::_checkerr(int errorcode)
{
  if (errorcode != 0)
    throw(PetscRuntimeError(errorcode));
}



inline
TiberPetscLinearSolver::~TiberPetscLinearSolver(void)
{
  this->clear();
}


inline
KSP TiberPetscLinearSolver::get_ksp(void)
{
  return _ksp;
}




#endif // _TIBERPETSCLINEARSOLVER_H_
