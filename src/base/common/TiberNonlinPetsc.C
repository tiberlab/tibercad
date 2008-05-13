// $Id$


#include "TiberNonlinPetsc.h"
#include "TiberPetscNonlinearSolver.h"

#include "InitFailedException.h"
#include "PetscDivergedError.h"


#include "equation_systems.h"
#include "mesh.h"


using namespace std;


TiberNonlinPetsc::TiberNonlinPetsc(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number)
{
  _solver = new TiberPetscNonlinearSolver<double>();

  if (_solver == NULL)
    throw InitFailedException("Cannot create linear solver object.");

}


TiberNonlinPetsc::~TiberNonlinPetsc(void)
{
  clear();

  delete _solver;
}




void
TiberNonlinPetsc::reinit(void)
{
  _solver->clear();
  _solver->init();

  Parent::reinit();
}



void
TiberNonlinPetsc::clear(void)
{
  _solver->clear();

  Parent::clear();
}



void
TiberNonlinPetsc::setup_pc(void)
{
  const string& pc = get_options()["pc_type"];
  if ((pc == "") || (pc == "ilu"))
    _solver->set_pc_type(PCILU);
  else if (pc == "lu")
    _solver->set_pc_type(PCLU);
  else if (pc == "jacobi")
    _solver->set_pc_type(PCJACOBI);
  else if (pc == "composite")
    _solver->set_pc_type(PCCOMPOSITE);
  else if (pc == "none")
    _solver->set_pc_type(PCNONE);
  else if (pc == "cholesky")
    _solver->set_pc_type(PCCHOLESKY);
  else
  {
    cerr << "PETSc nonlinear solver: unknown preconditioner \'"
      << pc << "\'. Falling back to \'ilu\'" << endl;
    _solver->set_pc_type(PCILU);
  }
}



void
TiberNonlinPetsc::setup_ksp(void)
{
  const string& ksp = get_options()["ksp_type"];
  if ((ksp == "") || (ksp == "bcgsl"))
    _solver->set_ksp_type(KSPBCGSL);
  else if (ksp == "bcgs")
    _solver->set_ksp_type(KSPBCGS);
  else if (ksp == "gmres")
    _solver->set_ksp_type(KSPGMRES);
  else if (ksp == "pconly")
    _solver->set_ksp_type(KSPPREONLY);
  else
  {
    cerr << "PETSc nonlinear solver: unknown Krylov method \'"
      << ksp << "\'. Falling back to \'bcgsl\'" << endl;
    _solver->set_ksp_type(KSPBCGSL);
  }
}



void
TiberNonlinPetsc::solve(void)
{

  assert(_assemble != NULL);

  _solver->matvec = _assemble;
  
  // Petsc uses l2 norms, but we specify tolerances per node
  // ('mean value' in the sense of eps = ||x|| / nn )
  double sqrt_nn = std::sqrt((double) get_mesh().n_nodes() * n_vars());

  _solver->set_ksp_options(_lin_tol, _lin_abs_tol, _lin_max_it);

  _solver->set_snes_options(_nonlin_rel_tol, _nonlin_abs_tol,
      _nonlin_step_tol, _nonlin_max_it);

  _solver->set_snes_ls_options(3, _max_step_size * sqrt_nn);

  KSPType ksp;
  switch (_solver_type)
  {
    case BICG:
      ksp = KSPBCGS;
      break;
      
    case GMRES:
      ksp = KSPGMRES;
      break;
      
    case BICGSTAB:
    default:
      ksp = KSPBCGSL;
      break;
  }
  _solver->set_ksp_type(ksp);
  
  PCType pc;
  switch (_preconditioner_type)
  {
    case IDENTITY_PRECOND:
      pc = PCNONE;
      break;

    case JACOBI_PRECOND:
      pc = PCJACOBI;
      break;

    case USER_PRECOND:
      pc = PCCOMPOSITE;
      break;

    case LU_PRECOND:
      pc = PCLU;
      break;

    case ILU_PRECOND:
    default:
      pc = PCILU;
      break;
  }
  _solver->set_pc_type(pc);


  bool failure = true;

  pair<unsigned int, double> result;

  try
  {
    result = _solver->solve(*matrix, *solution, *rhs, _lin_tol, _lin_max_it);

    failure = false;
  }
  catch (PetscDivergedError& e)
  {
    if (e.get_solver_type() == 1) cerr << "KSP ";
    else cerr << "SNES ";
    cerr << "diverged: " << e.get_reason() <<
      " at iteration " << e.get_iteration() <<
      " (fnorm = " << e.get_fnorm() << ")\n";

    //if (e.get_reason() == -5) retry = false;
    //if (e.get_reason() == -8) retry = false;
    //if (e.get_reason() == -6)
    //  solver_params.ls_type = 0;

    reinit();
    throw e;
  }
  catch (PetscRuntimeError& e)
  {
    std::cerr << "Petsc runtime error: " << e.get_reason() << std::endl;
    if (e.get_reason() == PETSC_ERR_MAT_LU_ZRPVT)
      cerr << " (Zero pivot during ILU.)\n";

    reinit();
    throw e;
  }


  _n_nonlin_iterations = result.first;
  _final_residual_norm = result.second;

  cout << "iterations: " << _n_nonlin_iterations <<
    ", residual = " << _final_residual_norm << endl;

  update();
}


