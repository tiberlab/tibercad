// $Id$

#include "TiberLinearSystem.h"
#include "TiberLinearSolver.h"


#include "equation_systems.h"
#include "linear_solver.h"

using namespace std;

TiberLinearSystem::TiberLinearSystem(EquationSystems& es,
    const string& name, const unsigned int number)
: TiberEqSystem(),
  LinearImplicitSystem(es, name, number)
{
  set_type(LINEAR);
}




TiberLinearSystem*
TiberLinearSystem::create(EquationSystems& es,
    const string& sysname, const ModelOptions& options)
{
  TiberLinearSystem* sys = NULL;
  sys = &(es.add_system<TiberLinearSystem>(sysname));

  assert(sys != NULL);
  sys->set_options(options);

  return sys;
}




void
TiberLinearSystem::user_initialization(void)
{
  TiberLinearSolver* lin_solver = TiberLinearSolver::create(get_options());
  linear_solver = AutoPtr<LinearSolver<Real> >(lin_solver);
}



void
TiberLinearSystem::solve(void)
{
  if (this->assemble_before_solve)
    this->assemble(); 

  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();
  
  TiberLinearSolver* lin_solver =
    static_cast<TiberLinearSolver*>(linear_solver.get());
  lin_solver->set_options(get_options());

  double lin_rel_tol = lin_solver->get_linear_rtol();
  int lin_max_it = lin_solver->get_linear_max_it();

  // Solve the linear system
  const std::pair<unsigned int, Real> rval;
  if (this->have_matrix("Preconditioner"))
    linear_solver->solve(*matrix, this->get_matrix("Preconditioner"),
        *solution, *rhs, lin_rel_tol, lin_max_it);
  else
    linear_solver->solve(*matrix, *solution, *rhs, lin_rel_tol, lin_max_it);


  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Update the system after the solve
  this->update();  
}
