// $Id$

#include "tiber_config.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#ifdef ENABLE_PARDISO
# include "PardisoLinearSolver.h"
#endif
#include "InitFailedException.h"
#include "ModelOptions.h"


namespace
{
  const double default_linear_rtol = 1e-6;
  const double default_linear_atol = 1e-50;
  const int default_linear_max_it = 500;
}




TiberLinearSolver::TiberLinearSolver(void)
  : _linear_rtol(default_linear_rtol),
    _linear_atol(default_linear_atol),
    _linear_max_it(default_linear_max_it)
{
}



TiberLinearSolver*
TiberLinearSolver::create(const std::string& type)
{
  ModelOptions opts;
  opts["linear_solver"] = type;
  
  return create(opts);
}
 

TiberLinearSolver*
TiberLinearSolver::create(const ModelOptions& options)
{
  TiberLinearSolver* solver = NULL;

  std::string type(options.get_option("linear_solver", "petsc"));

  if (type == "petsc")
    solver = new TiberPetscLinearSolver();
#ifdef ENABLE_PARDISO
  else if (type == "pardiso")
    solver = new PardisoLinearSolver();
#endif

  if (solver == NULL)
  {
    std::string msg = "TiberLinearSolver: no such solver '";
    msg += type + "'";
    throw InitFailedException(msg);
  }
#ifdef DEBUG
  std::cerr << "Created linear solver type " << type << std::endl;
#endif

  return solver;
}
         

void
TiberLinearSolver::set_options(const ModelOptions& options)
{
  _linear_rtol = options.get_option("lin_rel_tol", default_linear_rtol);
  _linear_atol = options.get_option("lin_abs_tol", default_linear_atol);
  _linear_max_it = options.get_option("lin_max_it", default_linear_max_it);

  _sim_name = options.get_option("name", "?");

  parse_options(options);
}
