// $Id$


#ifndef _TIBERPETSCLINEARSOLVER_H_
#define _TIBERPETSCLINEARSOLVER_H_

#include "TiberLinearSolver.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"

// Libmesh includes
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_macro.h"


EXTERN_C_FOR_PETSC_BEGIN
# include <petscversion.h>
# include <petscksp.h>
EXTERN_C_FOR_PETSC_END


//! The TiberCAD linear solver interface to PETSc
/*!
 * This class provides the TiberCAD interface to PETSc
 * iterative solvers.
 * It is derived from the libmesh LinearSolver class
 *
 */
class TiberPetscLinearSolver :  public TiberLinearSolver
{

  public:

    //!  Constructor. Initializes Petsc data structures
    TiberPetscLinearSolver(const libMesh::Parallel::Communicator &comm_in, const ModelOptions& options);


    //! Destructor.
    virtual ~TiberPetscLinearSolver(void);


    //! Release all memory and clear data structures.
    virtual void clear(void) override;


    //! Initialize data structures if not done so already.
    virtual void init(const char* name = NULL) override;


    //! Dummy implementation
    virtual std::pair<unsigned int, Real>
      solve(const libMesh::ShellMatrix<Number>&,
          const SparseMatrix<Number>&,
          NumericVector<Number>&,
          NumericVector<Number>&, double, unsigned int) override
      {
        return(std::make_pair(0, 0.0));
      };


    //! Dummy implementation
    virtual std::pair<unsigned int, Real>
      solve(const libMesh::ShellMatrix<Number>&,
          NumericVector<Number>&,
          NumericVector<Number>&, double, unsigned int) override 
      {
        return(std::make_pair(0, 0.0));
      };


    //! Return the convergence reason
    virtual libMesh::LinearConvergenceReason get_converged_reason() const override;


    //! Dummy implementation
    virtual void print_converged_reason(void) override {};


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
    virtual void do_parse_options(void) override;



    //! Solve the linear system
    virtual std::pair<unsigned int, Real>
      do_solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs) override;


    //! Setup the textual and graphical convergence monitors
    void setup_monitors(void);



  private:


    //! Krylov subspace context
    KSP _ksp{nullptr};

    //! The KSP type
    std::string _ksp_type{KSPBCGS};

    //! The PC type
    std::string _pc_type{PCILU};

    //! The graphical monitor
    PetscDrawLG _LG_monitor{nullptr};

    //! Do we want monitor?
    bool _monitor{false};

    //! Do we want X monitor?
    bool _xmonitor{false};

    //! Tells if xmonitor is already opened
    bool _xmonitor_open{false};

    //! The solver package
    std::string _solver_package{"petsc"};

    //! Check convergence
    std::pair<unsigned int, double> check_convergence(void);

    //! Set the sub PC
    void _set_sub_pc(PC pc, const std::string& pc_type);


};



//
// inline members
//



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
