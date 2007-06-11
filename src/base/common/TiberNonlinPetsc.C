// $Id$


#include "TiberNonlinPetsc.h"
#include "TiberPetscNonlinearSolver.h"

#include "InitFailedException.h"
#include "PetscDivergedError.h"
#include "SolveFailedException.h"


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
TiberNonlinPetsc::solve(void)
{

  assert(_assemble != NULL);

  _solver->matvec = _assemble;
  
  // Petsc uses l2 norms, but we specify tolerances per node
  // ('mean value' in the sense of eps = ||x|| / nn )
  double sqrt_nn = std::sqrt(get_mesh().n_nodes() * n_vars());

  _solver->set_ksp_options(_lin_tol, _lin_abs_tol, _lin_max_it);

  _solver->set_snes_options(_nonlin_rel_tol, _nonlin_abs_tol,
      _nonlin_step_tol, _nonlin_max_it);

  _solver->set_snes_ls_options(3, _max_step_size * sqrt_nn);

  KSPType ksp;
  switch (_solver_type)
  {
    case BICGSTAB:
      ksp = KSPBCGSL;
      break;
    case BICG:
      ksp = KSPBCGS;
      break;
    case GMRES:
      ksp = KSPGMRES;
      break;
  }
  _solver->set_ksp_type(ksp);
  
  PCType pc;
  switch (_preconditioner_type)
  {
    case ILU_PRECOND:
      pc = PCILU;
      break;
    case JACOBI_PRECOND:
      pc = PCJACOBI;
      break;
    case USER_PRECOND:
      pc = PCCOMPOSITE;
      break;
  }
  _solver->set_pc_type(pc);


  bool failure = true;
  string msg("DriftDiffusion: solve failed (");

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

    msg += e.what();
    msg += ")\n";

  }
  catch (PetscRuntimeError& e)
  {
    std::cerr << "Petsc runtime error: " << e.get_reason() << std::endl;
    if (e.get_reason() == PETSC_ERR_MAT_LU_ZRPVT)
    {
      // in the case of a zero pivot in (I)LU factorization
      // we try another preconditioner
      cerr << " (Zero pivot during ILU.)\n";
      try
      {
        _solver->set_pc_type(PCCOMPOSITE);
        result = _solver->solve(*matrix, *solution, *rhs, _lin_tol, _lin_max_it);
        failure = false;
      }
      catch (...)
      {
      }
    }
  }

  if (failure)
  {
    // we rebuild the equation system as could have been
    // 'damaged' by the crash
    reinit();

    throw SolveFailedException(msg);
  }

  _n_nonlin_iterations = result.first;
  _final_residual_norm = result.second;

  cerr << "iterations: " << _n_nonlin_iterations <<
    ", |r| = " << _final_residual_norm << endl;

  update();
}


