/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TiberPetscLinearSolver.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef TC_TIBERPETSCLINEARSOLVER_H
#define TC_TIBERPETSCLINEARSOLVER_H

#include "solver/TiberLinearSolver.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"

// Libmesh includes
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_macro.h"


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




#endif // TC_TIBERPETSCLINEARSOLVER_H
