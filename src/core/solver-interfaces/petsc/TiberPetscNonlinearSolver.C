// $Id$

#include "TiberPetscNonlinearSolver.h"
#include "TiberPetscUtils.h"
#include "ModelOptions.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"

// libmesh includes
#include "libmesh_common.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"

// C++ includes

#define DEBUG


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
      ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)))
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

#ifdef DEBUG
# if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR < 3))) 
    if (its == 0)
      std::cerr << "it " << its << ", fnorm = " << fnorm << "\n";
# endif
#endif

    if (fnorm != fnorm)
    {
#ifndef DEBUG
      std::cout << std::endl << std::flush;
#endif
      throw(KSPDivergedError(-8, its, fnorm));
    }

    // check for convergence
    //if ((reason < 0) && (reason != -3) && (reason != -4))
    if (reason < 0)
    //if ((reason < 0) && (reason != -3))
    {
#ifndef DEBUG
      std::cout << std::endl << std::flush;
#endif
      throw(KSPDivergedError(reason, its, fnorm));
    }

#ifdef DEBUG
# if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR < 3))) 
    if (its == 0)
      std::cerr << "it " << its << ", fnorm = " << fnorm << "\n";
# endif
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
    
    PetscVector<Number> X_global(x), R(r);
    PetscVector<Number> X_local(X_global.size());

    X_global.localize (X_local);
  
    if (solver->residual != NULL) solver->residual(X_local, R);
    if (solver->matvec   != NULL) solver->matvec(X_local, &R, NULL);

    R.close();
    
    return ierr;
  }


  
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __tiber_petsc_snes_jacobian(SNES snes, Vec x, Mat *jac, Mat *pc,
      MatStructure *msflag, void *ctx)
  {
    int ierr=0;

    assert (ctx != NULL);
    
    TiberPetscNonlinearSolver* solver =
      static_cast<TiberPetscNonlinearSolver*>(ctx);
    
    PetscMatrix<Number> PC(*pc);
    PetscMatrix<Number> Jac(*jac);
    PetscVector<Number> X_global(x);
    PetscVector<Number> X_local(X_global.size());

    X_global.localize(X_local);

    if (solver->jacobian != NULL) solver->jacobian(X_local, PC);
    if (solver->matvec   != NULL) solver->matvec(X_local, NULL, &PC);
    
    PC.close();
    Jac.close();
    
    *msflag = SAME_NONZERO_PATTERN;
    //*msflag = DIFFERENT_NONZERO_PATTERN;
    
    return ierr;
  }


  // This is the PETSc error handler
  PetscErrorCode
  __tiber_petsc_snes_error_handler(int line, const char* func,
      const char* file, const char* dir, PetscErrorCode n,
      int p, const char* mess, void* ctx) throw (PetscRuntimeError)
  {
    throw(PetscRuntimeError(n));
  }





  //
  // Convergence test
  // 
# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2))
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscInt it, PetscReal xnorm,
      PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

#ifdef DEBUG
    std::cerr << "iteration " << it << ": step = " << gnorm <<
      " residual = " << fnorm << std::endl;
#else
    std::cout << "." << std::flush;
#endif

    // this is a somewhat primitive check for divergence based on the step
    // norm (gnorm)
    if (it > 1)
    {
      TiberPetscNonlinearSolver* solver =
        static_cast<TiberPetscNonlinearSolver*>(ctx);
      if (gnorm > solver->get_divergence_tol() * solver->old_gnorm())
        throw(SNESDivergedError(-1, it, fnorm));

      solver->old_gnorm() = gnorm;
    }

    return SNESConverged_LS(snes, it, xnorm, gnorm, fnorm, reason, ctx);
    // TODO check gnorm
  }
#else
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscReal xnorm, PetscReal gnorm, 
      PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

#ifdef DEBUG
    std::cerr << "xnorm = " << xnorm << " gnorm = " << gnorm <<
      " fnorm = " << fnorm << "\n";
#else
    std::cout << "." << std::flush;
#endif

    return SNESConverged_LS(snes, xnorm, gnorm, fnorm, reason, ctx);
  }
#endif
    
} // end extern "C"







TiberPetscNonlinearSolver::TiberPetscNonlinearSolver(void)
  : _emergency_fnorm(1e-3),
    _ls_type(3),
    _ls_maxstep(1e5),
    _old_gnorm(1e96),
    _divergence_tol(2.0),
    _ksp_type(KSPBCGSL),
    _pc_type(PCILU)
{
}


void
TiberPetscNonlinearSolver::parse_options(const ModelOptions& options)
{
  _ls_type = TiberPetscUtils::extract_LSType(options);

  _ksp_type = TiberPetscUtils::extract_KSPType(options);

  _pc_type = TiberPetscUtils::extract_PCType(options);

  _ls_maxstep = options.get_option("max_step", _ls_maxstep);
  _divergence_tol = options.get_option("divergence_tol", _divergence_tol);
}



void TiberPetscNonlinearSolver::clear(void)
{
  if (this->initialized())
  {
    this->_is_initialized = false;

    int ierr=0;

    ierr = SNESDestroy(_snes);
    _checkerr(ierr);
  }
}



void TiberPetscNonlinearSolver::init(void)
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
    ierr = SNESCreate(libMesh::COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &_snes);
    _checkerr(ierr);

#else

    ierr = SNESCreate(libMesh::COMM_WORLD,&_snes);
    _checkerr(ierr);

#endif	     

#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR > 2))
    SNESMonitorSet(_snes, __tiber_petsc_snes_monitor, (void*) this, PETSC_NULL);
#else
    SNESSetMonitor(_snes, __tiber_petsc_snes_monitor, (void*) this, PETSC_NULL);
#endif

    PetscPushErrorHandler(__tiber_petsc_snes_error_handler, (void*) this);


    ierr = SNESSetType(_snes, SNESLS);
    _checkerr(ierr);


    SNESSetConvergenceTest(_snes, __tiber_snes_convergence_test, (void*) this);

    KSP ksp;
    SNESGetKSP(_snes, &ksp);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);

  }
}




std::pair<unsigned int, Real> 
TiberPetscNonlinearSolver::solve(SparseMatrix<double>&  jacobian,
    NumericVector<double>& solution,
    NumericVector<double>& residual)
{

  PetscMatrix<double>* jac = dynamic_cast<PetscMatrix<double>*>(&jacobian);
  PetscVector<double>* x   = dynamic_cast<PetscVector<double>*>(&solution);
  PetscVector<double>* r   = dynamic_cast<PetscVector<double>*>(&residual);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(jac != NULL); assert(jac->mat() != NULL);
  assert(x   != NULL); assert(x->vec()   != NULL);
  assert(r   != NULL); assert(r->vec()   != NULL);

  int ierr = 0;
  int n_iterations = 0;

  // for the case it was cleared before
  this->init();

  // setup the old_gnorm for divergence test
  _old_gnorm = 1e96;

  // set solver options
  SNESSetTolerances(_snes, get_nonlinear_atol(), get_nonlinear_rtol(),
      get_nonlinear_stol(), get_nonlinear_max_it(), get_linear_max_it());

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2))
  SNESSetMaxLinearSolveFailures(_snes, get_nonlinear_max_it());
#endif

  switch (_ls_type)
  {
    case 1:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESNoLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSet(_snes, SNESLineSearchNo, (void*) this);
#endif
      break;
    case 2:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESQuadraticLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSet(_snes, SNESLineSearchQuadratic, (void*) this);
#endif
      break;
    default:
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      ierr = SNESSetLineSearch(_snes, SNESCubicLineSearch, (void*) this);
#else
      ierr = SNESLineSearchSet(_snes, SNESLineSearchCubic, (void*) this);
#endif
      break;
  }

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  SNESSetLineSearchParams(_snes, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT);
  //SNESSetLineSearchParams(_snes, 0.00001, _ls_maxstep, PETSC_DEFAULT);
#else
  SNESLineSearchSetParams(_snes, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT);
#endif
  

  KSP ksp;
  SNESGetKSP(_snes, &ksp);

  ierr = KSPSetType(ksp, _ksp_type);
  _checkerr(ierr);

  KSPSetTolerances(ksp, get_linear_rtol(), get_linear_atol(), PETSC_DEFAULT,
      get_linear_max_it());
 
  PC pc;
  KSPGetPC(ksp, &pc);

  // get the type of preconditioner
  PCType pc_type = 0;
  PCGetType(pc, &pc_type);

  // - the very first time, there's no preconditioner yet
  // - if we changed the preconditioner, then we create it from scratch
  if ((pc_type == NULL) || (strcmp(_pc_type, pc_type) != 0))
  {
    ierr = PCSetType(pc, _pc_type);
    _checkerr(ierr);
    PCGetType(pc, &pc_type);

    // for composite type, do some extra stuff
    if (strcmp(pc_type, PCCOMPOSITE) == 0)
    {
      ierr = PCCompositeSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
      _checkerr(ierr);

      ierr = PCCompositeAddPC(pc, PCJACOBI);
      _checkerr(ierr);
      ierr = PCCompositeAddPC(pc, PCILU);
      _checkerr(ierr);

      PC sub_pc;
      ierr = PCCompositeGetPC(pc, 1, &sub_pc);
      _checkerr(ierr);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      PCILUSetZeroPivot(sub_pc, 1e-32);
      PCILUSetDamping(sub_pc, 1e-3);
#else
      PCFactorSetShiftNonzero(sub_pc, 1e-3);
      PCFactorSetZeroPivot(sub_pc, 1e-32);
      //PCILUReorderForNonzeroDiagonal(sub_pc, 1e-32);
#endif
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
  //_checkerr(ierr);


  // set functions
  SNESSetFunction (_snes, r->vec(),
      __tiber_petsc_snes_residual, (void*) this);

  SNESSetJacobian (_snes, jac->mat(), jac->mat(),
      __tiber_petsc_snes_jacobian, (void*) this);

  // solve the system
  // the syntax depends on the PETSc version
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

  ierr = SNESSolve(_snes, x->vec(), &n_iterations);
  _checkerr(ierr);

  // 2.2.x style	
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

  ierr = SNESSolve(_snes, x->vec());
  _checkerr(ierr);

  // 2.3.x & newer style	
#else

  ierr = SNESSolve(_snes, PETSC_NULL, x->vec());
  _checkerr(ierr);

#endif

  // check for convergence
  SNESConvergedReason reason;
  SNESGetConvergedReason(_snes, &reason);
  
  SNESGetIterationNumber(_snes, &n_iterations);

  double fnorm;
  SNESGetFunctionNorm(_snes, &fnorm);

  // reason < 0 means that the solver diverged. In this case we
  // throw an exception.
  if (reason < 0)
  { 
# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2))
    bool throw_ex = true;
    if (reason == -3)
    {
      KSPConvergedReason ksp_reason;
      KSPGetConvergedReason(ksp, &ksp_reason);
      if (ksp_reason == -3 || ksp_reason >= 0)
        throw_ex = false;
    }
    else if ((reason == -6) && (n_iterations == 1) && (fnorm < _emergency_fnorm))
      throw_ex = false;

    if (throw_ex)
      throw (SNESDivergedError(reason, n_iterations, fnorm));
#else
    if (reason != -3)
      if (!((reason == -6) && (n_iterations == 1) && (fnorm < _emergency_fnorm)))
        throw (SNESDivergedError(reason, n_iterations, fnorm));
#endif
  }

  std::cout << std::endl;

  // return the # of its. and the final residual norm.  Note that
  // n_iterations may be zero for PETSc versions 2.2.x and greater.
  return std::make_pair(n_iterations, fnorm);
}




