// $Id$


#include "TiberNonlinearSystem.h"

// the implemented systems/solvers
#include "TiberNonlinLS.h"
#include "TiberNonlinPetsc.h"

#include "InitFailedException.h"


#include "linear_solver.h"
#include "equation_systems.h"


using namespace std;


TiberNonlinearSystem::TiberNonlinearSystem(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number),
  _n_nonlin_iterations(0),
  _final_residual_norm(1e20),
  _last_step_size(1e20),
  _nonlin_step_tol(1e-3),
  _nonlin_rel_tol(1e-15),
  _nonlin_abs_tol(1e-15),
  _nonlin_max_it(20),
  _lin_tol(1e-6),
  _lin_max_it(500),
  _solver_type(BICGSTAB),
  _preconditioner_type(ILU_PRECOND)
{
}



TiberNonlinearSystem&
TiberNonlinearSystem::create_nonlinear_system(EquationSystems& es,
    const std::string& sysname, NonlinearSystemType type,
    const std::string& linear_solver)
{
  TiberNonlinearSystem* sys;
  
  switch (type)
  {
    case TIBER:
      sys = &(es.add_system<TiberNonlinLS>(sysname));
      break;
    default: // PETSc
      sys = &(es.add_system<TiberNonlinPetsc>(sysname));
      break;
  }

  sys->_linear_solver = linear_solver;
  return *sys;
}



TiberNonlinearSystem&
TiberNonlinearSystem::create_nonlinear_system(EquationSystems& es,
    const std::string& sysname, const std::string& type,
    const std::string& linear_solver)
{
  if (type == "tiber")
    return create_nonlinear_system(es, sysname, TIBER, linear_solver);
  else if (type == "petsc")
    return create_nonlinear_system(es, sysname, PETSC, linear_solver);
  else
  {
    std::string s = "TiberNonlinearSystem: unknown system type '" +
      type + "' for system " + sysname;
    throw InitFailedException(s);
  }
}
