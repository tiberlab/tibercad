// $Id$


#include "TiberPetscLinearSolver.h"

#include "libmesh_common.h"




void TiberPetscLinearSolver::clear(void)
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr = 0;
      
      ierr = KSPDestroy(_ksp);
      _checkerr(ierr);
	     
      // Mimic PETSc default solver and preconditioner
      this->_solver_type = GMRES;

      if (libMesh::n_processors() == 1)
	this->_preconditioner_type = ILU_PRECOND;
      else
	this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
}



void TiberPetscLinearSolver::init(void)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
  {
    this->_is_initialized = true;

    int ierr=0;


    // Create the linear solver context
    ierr = KSPCreate(libMesh::COMM_WORLD, &_ksp);
    _checkerr(ierr);

    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc);
    _checkerr(ierr);

    // We start with 0 for the correction
    //ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
    //_checkerr(ierr);

    // Set user-specified  solver and preconditioner types
    set_ksp_type();
    set_pc_type();

    // Set the options from user-input (for tests only)
    //ierr = KSPSetFromOptions (_ksp);
    //_checkerr(ierr);


    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
        PETSC_NULL,   // pointer to the array which holds the history
        PETSC_DECIDE, // size of the array holding the history
        PETSC_TRUE);  // Whether or not to reset the history for each solve. 
    _checkerr(ierr);
  }
}








std::pair<unsigned int, Real> 
TiberPetscLinearSolver::solve(SparseMatrix<Number>&  matrix_in,
			     SparseMatrix<Number>&  precond_in,
			     NumericVector<Number>& solution_in,
			     NumericVector<Number>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{
  this->init ();
  
  PetscMatrix<Number>* matrix   = dynamic_cast<PetscMatrix<Number>*>(&matrix_in);
  PetscMatrix<Number>* precond  = dynamic_cast<PetscMatrix<Number>*>(&precond_in);
  PetscVector<Number>* solution = dynamic_cast<PetscVector<Number>*>(&solution_in);
  PetscVector<Number>* rhs      = dynamic_cast<PetscVector<Number>*>(&rhs_in);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(matrix   != NULL);
  assert(precond  != NULL);
  assert(solution != NULL);
  assert(rhs      != NULL);
  
  int ierr = 0;
  int its = 0, max_its = static_cast<int>(m_its);
  PetscReal final_resid = 0.;

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();
  precond->close();
  solution->close();
  rhs->close();



  // we want PETSc >= 2.2.1 !
      
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			 SAME_NONZERO_PATTERN);
  _checkerr(ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances(_ksp, this->get_linear_rtol(), this->get_linear_atol(),
 			   PETSC_DEFAULT, this->get_linear_max_it());

  // Solve the linear system
  ierr = KSPSolve(_ksp, rhs->vec(), solution->vec());
  _checkerr(ierr);

	 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber(_ksp, &its);
	 
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm(_ksp, &final_resid);
  _checkerr(ierr);
	 

  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}



void TiberPetscLinearSolver::get_residual_history(std::vector<double>& hist)
{
  int ierr = 0;
  int its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  double* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  _checkerr(ierr);

  // Check for early return
  if (its == 0) return;
  
  // Create space to store the result
  hist.resize(its);

  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i)
    {
      hist[i] = *p;
      p++;
    }
}




Real TiberPetscLinearSolver::get_initial_residual(void)
{
  int ierr = 0;
  int its  = 0;

  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  double* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  _checkerr(ierr);

  // Check no residual history
  if (its == 0)
    {
      std::cerr << "No iterations have been performed, returning 0." << std::endl;
      return 0.;
    }

  // Otherwise, return the value pointed to by p.
  return *p;
}




void TiberPetscLinearSolver::set_ksp_type(void)
{
  int ierr = 0;
  
  switch (this->_solver_type)
    {

    case CG:
      ierr = KSPSetType (_ksp, (char*) KSPCG);
      _checkerr(ierr);
      break;

    case CR:
      ierr = KSPSetType (_ksp, (char*) KSPCR);
      _checkerr(ierr);
      break;

    case CGS:
      ierr = KSPSetType (_ksp, (char*) KSPCGS);
      _checkerr(ierr);
      break;

    case BICG:
      ierr = KSPSetType (_ksp, (char*) KSPBICG);
      _checkerr(ierr);
      break;

    case TCQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTCQMR);
      _checkerr(ierr);
      break;
 
    case TFQMR:
      ierr = KSPSetType (_ksp, (char*) KSPTFQMR);
      _checkerr(ierr);
      break;

    case LSQR:
      ierr = KSPSetType (_ksp, (char*) KSPLSQR);
      _checkerr(ierr);
      break;

    case BICGSTAB:
      ierr = KSPSetType (_ksp, (char*) KSPBCGS);
      _checkerr(ierr);
      break;

    case MINRES:
      ierr = KSPSetType (_ksp, (char*) KSPMINRES);
      _checkerr(ierr);
      break;

    case GMRES:
      ierr = KSPSetType (_ksp, (char*) KSPGMRES);
      _checkerr(ierr);
      break;

    case RICHARDSON:
      ierr = KSPSetType (_ksp, (char*) KSPRICHARDSON);
      _checkerr(ierr);
      break;

    case CHEBYSHEV: 
      ierr = KSPSetType (_ksp, (char*) KSPCHEBYCHEV);
      _checkerr(ierr);
      break;

    default:
      throw PetscRuntimeError(PETSC_ERR_ARG_UNKNOWN_TYPE);
    }
}








void TiberPetscLinearSolver::set_pc_type(void)
{
  int ierr = 0;
 
  switch (this->_preconditioner_type)
  {
    case IDENTITY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCNONE);
      _checkerr(ierr);
      break;

    case CHOLESKY_PRECOND:
      ierr = PCSetType (_pc, (char*) PCCHOLESKY);
      _checkerr(ierr);

    case ICC_PRECOND:
      ierr = PCSetType (_pc, (char*) PCICC);
      _checkerr(ierr);
      break;

    case ILU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCILU);
      _checkerr(ierr);
      break;

    case LU_PRECOND:
      ierr = PCSetType (_pc, (char*) PCLU);
      _checkerr(ierr);
      break;

    case ASM_PRECOND:
      ierr = PCSetType (_pc, (char*) PCASM);
      _checkerr(ierr);
      break;

    case JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCJACOBI);
      _checkerr(ierr);
      break;

    case BLOCK_JACOBI_PRECOND:
      ierr = PCSetType (_pc, (char*) PCBJACOBI);
      _checkerr(ierr);
      break;

    case SOR_PRECOND:
      ierr = PCSetType (_pc, (char*) PCSOR);
      _checkerr(ierr);
      break;

    case EISENSTAT_PRECOND:
      ierr = PCSetType (_pc, (char*) PCEISENSTAT);
      _checkerr(ierr);
      break;

    case USER_PRECOND:
      // abused for composite PC
      ierr = PCSetType (_pc, (char*) PCCOMPOSITE);
      _checkerr(ierr);
      ierr = PCCompositeSetType(_pc, PC_COMPOSITE_MULTIPLICATIVE);
      ierr = PCCompositeAddPC(_pc, PCJACOBI);
      ierr = PCCompositeAddPC(_pc, PCILU);

      PC sub_pc;
      ierr = PCCompositeGetPC(_pc, 1, &sub_pc);
      _checkerr(ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
      PCILUSetZeroPivot(sub_pc, 1e-54);
      PCLUSetZeroPivot(sub_pc, 1e-54);
      PCILUSetDamping(sub_pc, 1e-3);
#else
      PCFactorSetShiftNonzero(sub_pc, 1e-3);
      PCFactorSetZeroPivot(sub_pc, 1e-54);
      //PCILUReorderForNonzeroDiagonal(sub_pc, 1e-32);
#endif
      break;

    default:
      throw PetscRuntimeError(PETSC_ERR_ARG_UNKNOWN_TYPE);
  }

 
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  PCILUSetZeroPivot(_pc, 1e-54);
  PCLUSetZeroPivot(_pc, 1e-54);
#else
  PCFactorSetZeroPivot(_pc, 1e-54);
  //PCILUReorderForNonzeroDiagonal(pc, 1e-32);
#endif
}




// Explicit instantiation
//template class TiberPetscLinearSolver<Number>;
 


