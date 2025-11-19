// $Id$

#include "tibercad/solver/TiberLinearSystem.h"
#include "core/solver-interfaces/TiberLinearSolver.h"
#include "tibercad/base/InitFailedException.h"

#include "tibercad/base/TiberCad.h"


#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"

#include <cassert>

using namespace std;


TiberLinearSystem::TiberLinearSystem(libMesh::EquationSystems& es,
    const string& name, const unsigned int number)
: TiberEqSystem(),
  LinearImplicitSystem(es, name, number)
{
  set_type(LINEAR);
}




TiberLinearSystem*
TiberLinearSystem::create(libMesh::EquationSystems& es,
    const string& sysname, const ModelOptions& options)
{
  TiberLinearSystem* sys = nullptr;
  sys = &(es.add_system<TiberLinearSystem>(sysname));

  if (sys == nullptr)
    throw InitFailedException("Could not create linear system " + sysname);

  sys->set_options(options);

  return sys;
}




void
TiberLinearSystem::user_initialization(void)
{
  TiberLinearSolver* lin_solver = TiberLinearSolver::create(this->comm(), get_options());
  linear_solver = std::unique_ptr<libMesh::LinearSolver<Real> >(lin_solver);
}



void
TiberLinearSystem::solve(void)
{
  if (this->assemble_before_solve)
    this->assemble();

  // Get a reference to the EquationSystems
  const libMesh::EquationSystems& es =
    this->get_equation_systems();

  TiberLinearSolver* lin_solver =
    static_cast<TiberLinearSolver*>(linear_solver.get());
  lin_solver->set_options(get_options());

  double lin_rel_tol = lin_solver->get_linear_rtol();
  int lin_max_it = lin_solver->get_linear_max_it();

  // Solve the linear system
  const std::pair<unsigned int, Real> rval;
  if (this->have_matrix("Preconditioner"))
    lin_solver->solve(*matrix, this->get_matrix("Preconditioner"),
        *solution, *rhs, lin_rel_tol, lin_max_it);
  else
    lin_solver->solve(*matrix, *solution, *rhs, lin_rel_tol, lin_max_it);


  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;

  // Update the system after the solve
  this->update();

  lin_solver->clear();
}
