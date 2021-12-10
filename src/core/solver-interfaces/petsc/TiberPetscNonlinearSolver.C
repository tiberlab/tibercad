// $Id$

#include "TiberPetscNonlinearSolver.h"
#include "TiberPetscUtils.h"
#include "TiberCad.h"
#include "ModelOptions.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"

// libmesh includes
#include "libmesh_common.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "nonlinear_implicit_system.h"

#include "Messages.h"

#include <cstring>
#include <sstream>

// C++ includes
#include <cassert>
#include <cstring>


/*
 * These functions get called from the nonlinear solver of PETSc
 */
extern "C"
{


  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
# if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 2) \
      && (PETSC_VERSION_SUBMINOR == 0)) || \
      ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1))) \
	|| (PETSC_VERSION_MAJOR >= 3)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

  // this function is called by PETSc at the end of each nonlinear step
  PetscErrorCode
  __tiber_petsc_snes_monitor(SNES _snes, PetscInt its, PetscReal fnorm,
      void *ctx)
  {

    int ierr = 0;

    KSP ksp;
    SNESGetKSP(_snes, &ksp);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);

    if (fnorm != fnorm)
    {
      Messages::newline();
      throw(KSPDivergedError(-8, its, fnorm));
    }

    // check for convergence
    //if ((reason < 0) && (reason != -3) && (reason != -4))
    //if (reason < 0)
    if ((reason < 0) && (reason != -3))
    {
      Messages::newline();
      throw(KSPDivergedError(reason, its, fnorm));
    }

#if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR < 3)))
    if (its == 0)
    {
      std::ostringstream os;
      os << "it " << its << ", fnorm = " << fnorm;
      Messages::info(os.str());
    }
#endif

    TiberPetscNonlinearSolver* solver =
      static_cast<TiberPetscNonlinearSolver*>(ctx);
    solver->draw_point(its, fnorm);

    // we increase tolerance for linear solver at each nonlinear step
    double rtol, atol, stol;
    int maxit;
    KSPGetTolerances(ksp, &rtol, &atol, &stol, &maxit);
    KSPSetTolerances(ksp, rtol * rtol, atol / 100.0, stol, maxit);


    return ierr;
  }




  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __tiber_petsc_snes_residual(SNES, Vec x, Vec r, void *ctx)
  {
    int ierr=0;

    assert (x   != NULL);
    assert (r   != NULL);
    assert (ctx != NULL);

    TiberPetscNonlinearSolver* solver =
      static_cast<TiberPetscNonlinearSolver*>(ctx);

    libMesh::NonlinearImplicitSystem &sys = solver->system();

    libMesh::PetscVector<Number> X_global(x , solver->comm());
    libMesh::PetscVector<Number> R(r , solver->comm());
    libMesh::PetscVector<Number> X_local(solver->comm());
    X_local.init(X_global, true);

    // TODO see in libMesh for modifications

    X_global.localize (X_local);

    R.zero();

    if (solver->residual != NULL) solver->residual(X_local, R, sys);
    else if (solver->matvec   != NULL) solver->matvec(X_local, &R, NULL, sys);
    else if (solver->residual_and_jacobian_object != NULL)
      solver->residual_and_jacobian_object->residual_and_jacobian(X_local, &R, NULL, sys);

    R.close();
    X_global.close();

    return ierr;
  }



  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
# if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 5)
  __tiber_petsc_snes_jacobian(SNES, Vec x, Mat jac, Mat pc,
      void *ctx)
#else
  __tiber_petsc_snes_jacobian(SNES, Vec x, Mat *jac, Mat *pc,
      MatStructure *msflag, void *ctx)
#endif
  {
    int ierr=0;

    assert (ctx != NULL);

    TiberPetscNonlinearSolver* solver =
      static_cast<TiberPetscNonlinearSolver*>(ctx);

    libMesh::NonlinearImplicitSystem &sys = solver->system();

    libMesh::PetscMatrix<Number> PC(pc , solver->comm());
    libMesh::PetscMatrix<Number> Jac(jac , solver->comm());
    libMesh::PetscVector<Number> X_global(x , solver->comm());
    libMesh::PetscVector<Number> X_local(solver->comm());
    X_local.init(X_global, true);

    // Set the dof maps
    PC.attach_dof_map(sys.get_dof_map());
    Jac.attach_dof_map(sys.get_dof_map());

    X_global.localize(X_local);

    PC.zero();

    if (solver->jacobian != NULL) solver->jacobian(X_local, PC, sys);
    else if (solver->matvec   != NULL) solver->matvec(X_local, NULL, &PC, sys);
    else if (solver->residual_and_jacobian_object != NULL)
      solver->residual_and_jacobian_object->residual_and_jacobian(X_local, NULL, &PC, sys);

    PC.close();
    Jac.close();
    X_global.close();

#if PETSC_RELEASE_LESS_THAN(3,5,0)
    *msflag = SAME_NONZERO_PATTERN;
    //*msflag = DIFFERENT_NONZERO_PATTERN;
#endif

    return ierr;
  }


  // This is the PETSc error handler
  PetscErrorCode
  __tiber_petsc_snes_error_handler(MPI_Comm, int line, const char*,
# if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 5)
      const char* file, PetscErrorCode n,
#else
      const char* file, const char* dir, PetscErrorCode n,
#endif
      PetscErrorType, const char*, void*)
  {
#ifdef DEBUG
    std::ostringstream os;
    os << "PETSc error in " << file << ", line " << line;
    Messages::error(os.str());
#else
    (void*) file;
    (int) line;
#endif
    if ((n != 71) && (n != 81) && (n != 82))
      throw(PetscRuntimeError(n));
    else
      throw(SolverException("PETSc solver failed."));
  }





  //
  // Convergence test
  //
# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2)) || (PETSC_VERSION_MAJOR >= 3)
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscInt it, PetscReal xnorm,
      PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

    {
      std::ostringstream os;
      os << "iteration " << it << ": step = " << gnorm <<
        " residual = " << fnorm;
      Messages::info(os.str());
    }

    // this is a somewhat primitive check for divergence based on the step
    // norm (gnorm)
    if (it > 1)
    {
      TiberPetscNonlinearSolver* solver =
        static_cast<TiberPetscNonlinearSolver*>(ctx);
      if (gnorm > solver->get_divergence_tol() * solver->old_gnorm())
        throw(SNESDivergedError(-1, it, fnorm));

      solver->old_gnorm() = gnorm;

#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 4))
      return SNESConvergedDefault(snes, it, xnorm, gnorm, fnorm, reason, ctx);
#elif (PETSC_VERSION_MAJOR >= 3)
      return SNESDefaultConverged(snes, it, xnorm, gnorm, fnorm, reason, ctx);
#else
      return SNESConverged_LS(snes, it, xnorm, gnorm, fnorm, reason, ctx);
#endif
    }
    else
    {
      // this is a trick for situations where the solution is already found
      if (fnorm < 1e-12)
        *reason = SNES_CONVERGED_FNORM_ABS;
      return 0;
    }
    // TODO check gnorm
  }
#else
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscReal xnorm, PetscReal gnorm,
      PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

#ifdef DEBUG
    std::ostringstream os;
    os << "xnorm = " << xnorm << " gnorm = " << gnorm <<
      " fnorm = " << fnorm << "\n";
    Messages::info(os.str());
#else
    std::cout << "." << std::flush;
#endif

    return SNESConverged_LS(snes, xnorm, gnorm, fnorm, reason, ctx);
  }
#endif

} // end extern "C"






TiberPetscNonlinearSolver::TiberPetscNonlinearSolver(sys_type& s)
  : TiberNonlinearSolver(s),
    _emergency_fnorm(1e-3),
    _ls_type(3),
    _ls_maxstep(1e5),
    _old_gnorm(1e96),
    _divergence_tol(2.0),
    _ksp_type(KSPBCGSL),
    _pc_type(PCILU),
    _linear_rtol(1e-6),
    _linear_atol(1e-50),
    _linear_max_it(500)
{
}


void
TiberPetscNonlinearSolver::parse_options(const ModelOptions& options)
{
  _ls_type = TiberPetscUtils::extract_LSType(options);

  _ls_maxstep = options.get_option("max_step", _ls_maxstep);
  _divergence_tol = options.get_option("divergence_tolerance", _divergence_tol);

  ModelOptions::const_submodel_iterator it = options.submodels_begin("linear_solver");
  if (it != options.submodels_end("linear_solver"))
  {
    const ModelOptions& linoptions = it->second;
    _linear_rtol = linoptions.get_option("relative_tolerance", 1e-6);
    _linear_atol = linoptions.get_option("absolute_tolerance", 1e-50);
    _linear_max_it = linoptions.get_option("max_iterations", 500);

    _ksp_type = TiberPetscUtils::extract_KSPType(linoptions);
    _pc_type = TiberPetscUtils::extract_PCType(linoptions);
  }
}



void TiberPetscNonlinearSolver::clear(void)
{
  if (this->initialized())
  {
    this->_is_initialized = false;

    int ierr=0;

    ierr = SNESDestroy(&_snes);
    TiberPetscUtils::checkerr(ierr);
  }
}



void TiberPetscNonlinearSolver::init(const char*)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
  {
    this->_is_initialized = true;

    int ierr = 0;

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1) && \
    (PETSC_VERSION_SUBMINOR <= 1))

    // At least until Petsc 2.1.1, the SNESCreate had a different
    // calling syntax.
    // The second argument was of type SNESProblemType, and could have a
    // value of either SNES_NONLINEAR_EQUATIONS or
    // SNES_UNCONSTRAINED_MINIMIZATION.
    ierr = SNESCreate(comm().get(), SNES_NONLINEAR_EQUATIONS, &_snes);
    TiberPetscUtils::checkerr(ierr);

#else
    ierr = SNESCreate(comm().get(),&_snes);
    TiberPetscUtils::checkerr(ierr);

#endif

#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR > 2)) || (PETSC_VERSION_MAJOR >= 3)
    SNESMonitorSet(_snes, __tiber_petsc_snes_monitor, (void*) this, PETSC_NULL);
#else
    SNESSetMonitor(_snes, __tiber_petsc_snes_monitor, (void*) this, PETSC_NULL);
#endif

    PetscPushErrorHandler(__tiber_petsc_snes_error_handler, (void*) this);


    ierr = SNESSetType(_snes, SNESNEWTONLS);
    TiberPetscUtils::checkerr(ierr);


#if (PETSC_VERSION_MAJOR >= 3)
    SNESSetConvergenceTest(_snes, __tiber_snes_convergence_test,
        (void*) this, NULL);
#else
    SNESSetConvergenceTest(_snes, __tiber_snes_convergence_test,
        (void*) this);
#endif

    KSP ksp;
    SNESGetKSP(_snes, &ksp);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
  }
}




std::pair<unsigned int, Real>
TiberPetscNonlinearSolver::solve(SparseMatrix<double>&  jacobian,
    NumericVector<double>& solution,
    NumericVector<double>& residual)
{

  libMesh::PetscMatrix<double>* jac = dynamic_cast<libMesh::PetscMatrix<double>*>(&jacobian);
  libMesh::PetscVector<double>* x   = dynamic_cast<libMesh::PetscVector<double>*>(&solution);
  libMesh::PetscVector<double>* r   = dynamic_cast<libMesh::PetscVector<double>*>(&residual);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(jac != NULL); assert(jac->mat() != NULL);
  assert(x   != NULL); assert(x->vec()   != NULL);
  assert(r   != NULL); assert(r->vec()   != NULL);

  //x->close();
  //r->close();
  //jac->close();

  int ierr = 0;
  int n_iterations = 0;

  // for the case it was cleared before
  this->init();

  // setup the old_gnorm for divergence test
  _old_gnorm = 1e96;

  // set solver options
  SNESSetTolerances(_snes, get_nonlinear_atol(), get_nonlinear_rtol(),
      get_nonlinear_stol(), get_nonlinear_max_it(), _linear_max_it);

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2)) || (PETSC_VERSION_MAJOR >= 3)
  SNESSetMaxLinearSolveFailures(_snes, get_nonlinear_max_it());
#endif

  SNESLineSearch ls;
  SNESGetLineSearch(_snes, &ls);

  switch (_ls_type)
  {
    case 1:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESNoLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSetOrder(ls, SNES_LINESEARCH_ORDER_LINEAR);
#endif
      break;
    case 2:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESQuadraticLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSetOrder(ls, SNES_LINESEARCH_ORDER_QUADRATIC);
#endif
      break;
    default:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESCubicLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSetOrder(ls, SNES_LINESEARCH_ORDER_CUBIC);
#endif
      break;
  }

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  SNESSetLineSearchParams(_snes, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT);
  //SNESSetLineSearchParams(_snes, 0.00001, _ls_maxstep, PETSC_DEFAULT);
#elif (PETSC_VERSION_MAJOR >= 3)
  SNESLineSearchSetTolerances(ls, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
#else
  SNESLineSearchSetParams(_snes, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT);
#endif

  KSP ksp;
  SNESGetKSP(_snes, &ksp);

  ierr = KSPSetType(ksp, _ksp_type.c_str());
  TiberPetscUtils::checkerr(ierr);

  KSPSetTolerances(ksp, _linear_rtol, _linear_atol, PETSC_DEFAULT,
      _linear_max_it);

  PC pc;
  KSPGetPC(ksp, &pc);

  // get the type of preconditioner
  PCType pc_type = 0;
  PCGetType(pc, &pc_type);

  // - the very first time, there's no preconditioner yet
  // - if we changed the preconditioner, then we create it from scratch
  if ((pc_type == NULL) || (_pc_type.compare(pc_type) != 0))
  //if (0)
  {
    ierr = PCSetType(pc, _pc_type.c_str());
    TiberPetscUtils::checkerr(ierr);
    PCGetType(pc, &pc_type);

    // for composite type, do some extra stuff
    if (strcmp(pc_type, PCCOMPOSITE) == 0)
    {
      /*
      ierr = PCCompositeSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
      TiberPetscUtils::checkerr(ierr);

      ierr = PCCompositeAddPC(pc, (char*) PCJACOBI);
      TiberPetscUtils::checkerr(ierr);
      ierr = PCCompositeAddPC(pc, (char*) PCILU);
      TiberPetscUtils::checkerr(ierr);

      PC sub_pc;
      ierr = PCCompositeGetPC(pc, 1, &sub_pc);
      TiberPetscUtils::checkerr(ierr);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      PCILUSetZeroPivot(sub_pc, 1e-32);
      PCILUSetDamping(sub_pc, 1e-3);
#elif (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 1)
    PCFactorSetShiftType(sub_pc, MAT_SHIFT_NONZERO);
    PCFactorSetShiftAmount(sub_pc, 1e-3);
    PCFactorSetZeroPivot(sub_pc, 1e-32);
    //PCILUReorderForNonzeroDiagonal(sub_pc, 1e-32);
#else
      PCFactorSetShiftNonzero(sub_pc, 1e-3);
      PCFactorSetZeroPivot(sub_pc, 1e-32);
      //PCILUReorderForNonzeroDiagonal(sub_pc, 1e-32);
#endif
    */
    }
  }

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  PCILUSetZeroPivot(pc, 1e-32);
#else
  PCFactorSetZeroPivot(pc, 1e-32);
  //PCILUReorderForNonzeroDiagonal(pc, 1e-32);
#endif

  // to override options from command line
  // only for tests
  //ierr = SNESSetFromOptions(_snes);
  //TiberPetscUtils::checkerr(ierr);

  // set functions
  SNESSetFunction (_snes, r->vec(),
      __tiber_petsc_snes_residual, (void*) this);

  SNESSetJacobian (_snes, jac->mat(), jac->mat(),
      __tiber_petsc_snes_jacobian, (void*) this);

  // solve the system
  // the syntax depends on the PETSc version
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

  ierr = SNESSolve(_snes, x->vec(), &n_iterations);
  TiberPetscUtils::checkerr(ierr);

  // 2.2.x style
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

  ierr = SNESSolve(_snes, x->vec());
  TiberPetscUtils::checkerr(ierr);

  // 2.3.x & newer style
#else

  ierr = SNESSolve(_snes, PETSC_NULL, x->vec());
  TiberPetscUtils::checkerr(ierr);

#endif

  // check for convergence
  SNESConvergedReason reason;
  SNESGetConvergedReason(_snes, &reason);

  SNESGetIterationNumber(_snes, &n_iterations);

  double fnorm;
# if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 5))
  fnorm = r->l2_norm();
#else
  SNESGetFunctionNorm(_snes, &fnorm);
#endif

  // reason < 0 means that the solver diverged. In this case we
  // throw an exception.
  if (reason < 0)
  {
# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2)) || (PETSC_VERSION_MAJOR >= 3)
    bool throw_ex = true;
    if (reason == -3)
    {
      KSPConvergedReason ksp_reason;
      KSPGetConvergedReason(ksp, &ksp_reason);
      if (ksp_reason == -3 || ksp_reason >= 0)
        throw_ex = false;
    }
# if (PETSC_VERSION_MAJOR >= 3)
    else if ((reason == -6) && (n_iterations == 0) && (fnorm < _emergency_fnorm))
# else
    else if ((reason == -6) && (n_iterations == 1) && (fnorm < _emergency_fnorm))
# endif
      throw_ex = false;

    if (throw_ex)
      throw (SNESDivergedError(reason, n_iterations, fnorm));
#else
    if (reason != -3)
      if (!((reason == -6) && (n_iterations == 1) && (fnorm < _emergency_fnorm)))
        throw (SNESDivergedError(reason, n_iterations, fnorm));
#endif
  }

  // return the # of its. and the final residual norm.  Note that
  // n_iterations may be zero for PETSc versions 2.2.x and greater.
  return std::make_pair(n_iterations, fnorm);
}




