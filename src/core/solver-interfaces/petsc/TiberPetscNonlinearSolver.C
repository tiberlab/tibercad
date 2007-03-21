// $Id$


// libmesh includes
#include "libmesh_common.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"

// C++ includes

// Local Includes
#include "TiberPetscNonlinearSolver.h"



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
      void *) throw (PetscRuntimeError, KSPDivergedError)
  {

    int ierr = 0;

    KSP ksp;
    SNESGetKSP(_snes, &ksp);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);

    if (fnorm != fnorm)
      throw(KSPDivergedError(-8, its, fnorm));

    // check for convergence
    //if ((reason < 0) && (reason != -3) && (reason != -4))
    if ((reason < 0) && (reason != -3))
      throw(KSPDivergedError(reason, its, fnorm));

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
    
    TiberPetscNonlinearSolver<Number>* solver =
      static_cast<TiberPetscNonlinearSolver<Number>*>(ctx);
    
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
  __tiber_petsc_snes_jacobian(SNES, Vec x, Mat *jac, Mat *pc,
      MatStructure *msflag, void *ctx)
  {
    int ierr=0;

    assert (ctx != NULL);
    
    TiberPetscNonlinearSolver<Number>* solver =
      static_cast<TiberPetscNonlinearSolver<Number>*>(ctx);
    
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


  /*
  // Checks the update vector before updating
  PetscErrorCode
  __tiber_snes_update(SNES snes, PetscInt)
  {
    int ierr = 0;

    Vec x;

    ierr = SNESGetSolutionUpdate(snes, &x);
    if (ierr != 0) throw(PetscRuntimeError(ierr));
    
    PetscScalar* vv;
    ierr = VecGetArray(x, &vv);
    if (ierr != 0) throw(PetscRuntimeError(ierr));

    PetscInt n;
    ierr = VecGetLocalSize(x, &n);

    for (int i = 0; i < n; i++)
    {
      //if (vv[i] > 0.1)
      //{
        //vv[i] = 0.0001;
      //}
    }
    
    ierr = VecRestoreArray(x, &vv);
    if (ierr != 0) throw(PetscRuntimeError(ierr));
    
    return ierr;
  }
  */


  // Convergence test
# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2))
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscInt it, PetscReal xnorm,
      PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

    std::cerr << "iteration " << it << ": xnorm = " << xnorm <<
      " gnorm = " << gnorm << " fnorm = " << fnorm << "\n";

    return SNESConverged_LS(snes, it, xnorm, gnorm, fnorm, reason, ctx);
  }
#else
  PetscErrorCode
  __tiber_snes_convergence_test(SNES snes, PetscReal xnorm, PetscReal gnorm, 
      PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr = 0;

    std::cerr << "xnorm = " << xnorm << " gnorm = " << gnorm <<
      " fnorm = " << fnorm << "\n";

    return SNESConverged_LS(snes, xnorm, gnorm, fnorm, reason, ctx);
  }
#endif
    
} // end extern "C"







template <typename T>
TiberPetscNonlinearSolver<T>::TiberPetscNonlinearSolver(void)
  throw (PetscRuntimeError) 
  : _nonlinear_rtol(1e-15),
    _nonlinear_atol(1e-18),
    _nonlinear_stol(1e-6),
    _linear_rtol(1e-6),
    _linear_atol(1e-15),
    _nonlinear_max_it(20),
    _linear_max_it(500),
    _ls_type(3),
    _ls_maxstep(1e3),
    _ksp_type(KSPBCGSL),
    _pc_type(PCILU)
{
}


template <typename T>
void
TiberPetscNonlinearSolver<T>::set_snes_options(double rtol, double atol,
    double stol, unsigned int max_it)
{
  _nonlinear_rtol = rtol;
  _nonlinear_atol = atol;
  _nonlinear_stol = stol;
  _nonlinear_max_it = max_it;
}

template <typename T>
void
TiberPetscNonlinearSolver<T>::set_snes_ls_options(int ls_type, 
    double ls_maxstep)
{
  _ls_type = ls_type;
  _ls_maxstep = ls_maxstep;
}





template <typename T>
void TiberPetscNonlinearSolver<T>::clear(void) throw (PetscRuntimeError)
{
  if (this->initialized())
  {
    this->_is_initialized = false;

    int ierr=0;

    ierr = SNESDestroy(_snes);
    _checkerr(ierr);
  }
}



template <typename T>
void TiberPetscNonlinearSolver<T>::init(void) throw (PetscRuntimeError)
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

    SNESSetMonitor(_snes, __tiber_petsc_snes_monitor, (void*) this, PETSC_NULL);
    PetscPushErrorHandler(__tiber_petsc_snes_error_handler, (void*) this);

    //SNESSetUpdate(_snes, __tiber_snes_update);

    ierr = SNESSetType(_snes, SNESLS);
    _checkerr(ierr);


    SNESSetConvergenceTest(_snes, __tiber_snes_convergence_test, (void*) this);

    KSP ksp;
    SNESGetKSP(_snes, &ksp);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);

  }
}


template <typename T>
std::pair<unsigned int, Real> 
TiberPetscNonlinearSolver<T>::solve(SparseMatrix<T>&  jacobian,
    NumericVector<T>& solution,
    NumericVector<T>& residual,
    const double rtol,
    const unsigned int iter)
  throw (PetscRuntimeError, KSPDivergedError, SNESDivergedError)
{

  PetscMatrix<T>* jac = dynamic_cast<PetscMatrix<T>*>(&jacobian);
  PetscVector<T>* x   = dynamic_cast<PetscVector<T>*>(&solution);
  PetscVector<T>* r   = dynamic_cast<PetscVector<T>*>(&residual);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(jac != NULL); assert(jac->mat() != NULL);
  assert(x   != NULL); assert(x->vec()   != NULL);
  assert(r   != NULL); assert(r->vec()   != NULL);

  int ierr = 0;
  int n_iterations = 0;

  // for the case it was cleared before
  this->init();
  
  // set solver options
  SNESSetTolerances(_snes, _nonlinear_atol, _nonlinear_rtol,
      _nonlinear_stol, _nonlinear_max_it, _linear_max_it);

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && \
    (PETSC_VERSION_SUBMINOR >= 2))
  SNESSetMaxLinearSolveFailures(_snes, _nonlinear_max_it);
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
#else
  SNESLineSearchSetParams(_snes, PETSC_DEFAULT, _ls_maxstep, PETSC_DEFAULT);
#endif
  

  KSP ksp;
  SNESGetKSP(_snes, &ksp);

  ierr = KSPSetType(ksp, _ksp_type);
  _checkerr(ierr);
  
  KSPSetTolerances(ksp, _linear_rtol, _linear_atol, PETSC_DEFAULT,
      _linear_max_it);
 
  PC pc;
  KSPGetPC(ksp, &pc);

  // get the type of preconditioner
  PCType pc_type;
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
    if (throw_ex)
      throw (SNESDivergedError(reason, n_iterations, fnorm));
#else
    throw (SNESDivergedError(reason, n_iterations, fnorm));
#endif
  }

  // return the # of its. and the final residual norm.  Note that
  // n_iterations may be zero for PETSc versions 2.2.x and greater.
  return std::make_pair(n_iterations, fnorm);
}




// Explicit instantiation
template class TiberPetscNonlinearSolver<Real>;


