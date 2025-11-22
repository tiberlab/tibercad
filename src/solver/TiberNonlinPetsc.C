// $Id$


#include "TiberNonlinPetsc.h"

#include "solver/petsc/TiberPetscNonlinearSolver.h"

#include "tibercad/base/InitFailedException.h"
#include "solver/petsc/PetscDivergedError.h"
#include "tibercad/base/SolveFailedException.h"


#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh.h"

#include "tibercad/io/Messages.h"
#include <cassert>


using namespace std;


TiberNonlinPetsc::TiberNonlinPetsc(libMesh::EquationSystems& es,
    const std::string& name, const unsigned int number)
  : Parent(es, name, number)
{
  _solver = new TiberPetscNonlinearSolver(*this);

  if (_solver == NULL)
    throw InitFailedException("Cannot create nonlinear solver object.");

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
  _solver->set_xmonitor(NULL);

  Parent::reinit();
}



void
TiberNonlinPetsc::clear(void)
{
  _solver->clear();

  Parent::clear();
}




void
TiberNonlinPetsc::do_solve(void)
{

  if (_assemble != NULL)
    _solver->matvec = _assemble;
  else if (_solver->residual_and_jacobian_object == NULL)
  {
    throw SolveFailedException("No assembly routines set for nonlinear solver");
  }

  // setup the max line search step
  double sqrt_nn = std::sqrt((double) get_mesh().n_nodes() * n_vars());
  double ls_max_step = get_options().get_option("ls_max_step", 10.0);
  get_options().set_option("max_step", ls_max_step * sqrt_nn);

  _solver->set_options(get_options());
  
  bool failure = true;

  pair<unsigned int, double> result;

  try
  {
    _solver->set_xmonitor(get_xmonitor());
    result = _solver->solve(*matrix, *solution, *rhs);
    _solver->set_xmonitor(NULL);

    failure = false;
  }
  catch (PetscDivergedError& e)
  {
    //if (e.get_solver_type() == 1) cerr << "KSP ";
    //else cerr << "SNES ";
    //cerr << "diverged: " << e.get_reason() <<
    //  " at iteration " << e.get_iteration() <<
    //  " (fnorm = " << e.get_fnorm() << ")\n";

    //if (e.get_reason() == -5) retry = false;
    //if (e.get_reason() == -8) retry = false;
    //if (e.get_reason() == -6)
    //  solver_params.ls_type = 0;

    reinit();
    throw e;
  }
  catch (PetscRuntimeError& e)
  {
    //std::cerr << "Petsc runtime error: " << e.get_reason() << std::endl;
    //if (e.get_reason() == PETSC_ERR_MAT_LU_ZRPVT)
    //  cerr << " (Zero pivot during ILU.)\n";

    reinit();
    throw e;
  }


  _n_nonlin_iterations = result.first;
  _final_residual_norm = result.second;

  ostringstream os;
  os << "iterations: " << _n_nonlin_iterations <<
    ", final residual = " << _final_residual_norm;
  Messages::newline();
  Messages::info(os.str());

  update();
}


