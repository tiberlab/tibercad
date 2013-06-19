// $Id$


#include "TiberPetscLinearSolver.h"
#include "TiberPetscUtils.h"
#include "ModelOptions.h"

#include "libmesh_common.h"

#include <cassert>
#include <cstring>



TiberPetscLinearSolver::TiberPetscLinearSolver(const ModelOptions& options)
  : TiberLinearSolver(options),
    _ksp(NULL),
    _ksp_type(KSPBCGS),
    _pc_type(PCILU),
    _monitor(false),
    _xmonitor(true),
    _xmonitor_open(false)
{
//  if (libMesh::n_processors() == 1)
//    this->_preconditioner_type = ILU_PRECOND;
//  else
//    this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
}




void TiberPetscLinearSolver::clear(void)
{
  if (this->initialized())
    {
      _is_initialized = false;

      int ierr = 0;
      
      ierr = KSPDestroy(_ksp);
      TiberPetscUtils::checkerr(ierr);
	     
      // Mimic PETSc default solver and preconditioner
      //_solver_type = GMRES;

      //if (libMesh::n_processors() == 1)
      //  this->_preconditioner_type = ILU_PRECOND;
      //else
      //  this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
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
    TiberPetscUtils::checkerr(ierr);

    // Create the preconditioner context
    //ierr = KSPGetPC(_ksp, &_pc);
    //TiberPetscUtils::checkerr(ierr);

    // We start with 0 for the correction
    //ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
    //TiberPetscUtils::checkerr(ierr);

    // Set the options from user-input (for tests only)
    //ierr = KSPSetFromOptions (_ksp);
    //TiberPetscUtils::checkerr(ierr);


    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
        PETSC_NULL,   // pointer to the array which holds the history
        PETSC_DECIDE, // size of the array holding the history
        PETSC_TRUE);  // Whether or not to reset the history for each solve. 
    TiberPetscUtils::checkerr(ierr);
  }
}








std::pair<unsigned int, double> 
TiberPetscLinearSolver::do_solve(SparseMatrix<Number>&  matrix_in,
			     SparseMatrix<Number>&  precond_in,
			     NumericVector<Number>& solution_in,
			     NumericVector<Number>& rhs_in)
{
  this->init();
  
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

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();
  precond->close();
  solution->close();
  rhs->close();

  std::string ksp_type(_ksp_type);
  if ((_pc_type == PCLU) && (matrix == precond))
    ksp_type = KSPPREONLY;

  // Set user-specified solver and preconditioner types
  ierr = KSPSetType(_ksp, ksp_type.c_str());
  TiberPetscUtils::checkerr(ierr);

  PC pc;
  ierr = KSPGetPC(_ksp, &pc);
  TiberPetscUtils::checkerr(ierr);

  ierr = PCSetType(pc, _pc_type.c_str());

  if (_solver_package != "") {
    PCFactorSetMatSolverPackage(pc, _solver_package.c_str());
  }

  TiberPetscUtils::checkerr(ierr);
  // for composite type, do some extra stuff
  if (_pc_type.compare(PCCOMPOSITE) == 0)
  {
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
#else
    PCFactorSetShiftNonzero(sub_pc, 1e-3);
    PCFactorSetZeroPivot(sub_pc, 1e-32);
    //PCILUReorderForNonzeroDiagonal(sub_pc, 1e-32);
#endif
  }
 
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  PCILUSetZeroPivot(pc, 1e-54);
  PCLUSetZeroPivot(pc, 1e-54);
#else
  PCFactorSetZeroPivot(pc, 1e-54);
  //PCILUReorderForNonzeroDiagonal(pc, 1e-32);
#endif

  setup_monitors();

  // we want PETSc >= 2.2.1 !
      
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			 SAME_NONZERO_PATTERN);
  TiberPetscUtils::checkerr(ierr);

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances(_ksp, get_linear_rtol(), get_linear_atol(),
 			   PETSC_DEFAULT, get_linear_max_it());

  // Solve the linear system
  ierr = KSPSolve(_ksp, rhs->vec(), solution->vec());
  //KSPConvergedReason reason;
  //KSPGetConvergedReason(_ksp, &reason);
  //std::cerr << "KSP convergence reason: " << reason << std::endl;
  TiberPetscUtils::checkerr(ierr);

  return check_convergence();

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
  TiberPetscUtils::checkerr(ierr);

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




double TiberPetscLinearSolver::get_initial_residual(void)
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
  TiberPetscUtils::checkerr(ierr);

  // Check no residual history
  if (its == 0)
    {
      std::cerr << "No iterations have been performed, returning 0." << std::endl;
      return 0.;
    }

  // Otherwise, return the value pointed to by p.
  return *p;
}




std::pair<unsigned int, double>
TiberPetscLinearSolver::check_convergence(void)
{
	
  int its = 0;
  double fnorm = 0.0;
	 
  // Get the number of iterations required for convergence
  KSPGetIterationNumber(_ksp, &its);
	 
  // Get the norm of the final residual to return to the user.
  KSPGetResidualNorm(_ksp, &fnorm);

  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);

  if (reason <= 0)
  {
    bool throw_ex = true;

    if ((reason == -3) && (fnorm < 1e-9))
      throw_ex = false;

    if (throw_ex)
      throw(KSPDivergedError(reason, its, fnorm));
  }

  return std::pair<unsigned int, double>(its, fnorm);
}


void
TiberPetscLinearSolver::setup_monitors(void)
{
  if (_ksp != NULL)
  {
    int ierr = 0;
    if (_monitor)
    {
#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
    && (PETSC_VERSION_SUBMINOR > 2)) || (PETSC_VERSION_MAJOR >= 3)
      ierr = KSPMonitorSet(_ksp, KSPMonitorDefault, PETSC_NULL, 0);
#else
      ierr = KSPSetMonitor(_ksp, KSPDefaultMonitor, PETSC_NULL, 0);
#endif
    }
    else
      ierr = KSPMonitorCancel(_ksp);
      // does not work:
      //ierr = KSPMonitorSet(_ksp, PETSC_NULL, PETSC_NULL, 0);

    TiberPetscUtils::checkerr(ierr);


    if (_xmonitor)
    {
      if (!_xmonitor_open)
      {
        std::string sim_name(get_simulation_name());
        sim_name += ": Linear solver convergence monitor";
#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3)	\
    && (PETSC_VERSION_SUBMINOR > 2)) || (PETSC_VERSION_MAJOR >= 3)
        ierr = KSPMonitorLGCreate(NULL, sim_name.c_str(),0,0,400,400, &_LG_monitor);
#else
        ierr = KSPLGMonitorCreate(NULL, sim_name.c_str(),0,0,400,400, &_LG_monitor);
#endif
        TiberPetscUtils::checkerr(ierr);
    
#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3)	\
    && (PETSC_VERSION_SUBMINOR > 2)) || (PETSC_VERSION_MAJOR >= 3)
        ierr = KSPMonitorSet(_ksp, KSPMonitorLG, _LG_monitor, 0);
#else
        ierr = KSPSetMonitor(_ksp, KSPLGMonitor, _LG_monitor, 0);
#endif
        TiberPetscUtils::checkerr(ierr);
      
        _xmonitor_open = true;
      }
    }
    else if (_xmonitor_open)
    {
      // close xmonitor
    }
  }
}



void
TiberPetscLinearSolver::do_parse_options(void)
{
  _ksp_type = TiberPetscUtils::extract_KSPType(get_options());

  _pc_type = TiberPetscUtils::extract_PCType(get_options());

  _solver_package = get_option("solver_package", "");

  _monitor = get_option("monitor", false);
  _xmonitor = get_option("xmonitor", false);
}



