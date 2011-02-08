// $Id$

#include "tiber_config.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#ifdef ENABLE_PARDISO
# include "PardisoLinearSolver.h"
#endif
#include "InitFailedException.h"
#include "ModelOptions.h"
#include "Messages.h"


namespace
{
  const double default_linear_rtol = 1e-6;
  const double default_linear_atol = 1e-50;
  const int default_linear_max_it = 1500;
}




TiberLinearSolver::TiberLinearSolver(const ModelOptions& options)
  : TiberModelObject(options),
    _linear_rtol(default_linear_rtol),
    _linear_atol(default_linear_atol),
    _linear_max_it(default_linear_max_it)
{
}




TiberLinearSolver*
TiberLinearSolver::create(const ModelOptions& options)
{
  TiberLinearSolver* solver = NULL;

  std::string type(options.get_name());
  if (type.empty())
    type = "petsc";
  type = options.get_option("type", type);

  if (type == "petsc")
    solver = new TiberPetscLinearSolver(options);
#ifdef ENABLE_PARDISO
  else if (type == "pardiso")
    solver = new PardisoLinearSolver(options);
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

  // dummy read
  solver->get_option("type", "");
  solver->get_option("simulation", "");

  return solver;
}




void
TiberLinearSolver::parse_options(void)
{
  _linear_rtol = get_option("relative_tolerance", default_linear_rtol);
  _linear_atol = get_option("absolute_tolerance", default_linear_atol);
  _linear_max_it = get_option("max_iterations", default_linear_max_it);

  do_parse_options();

  get_options().check_unused();
}



const std::string&
TiberLinearSolver::get_simulation_name(void) const
{
  return get_option("simulation", "");
}


std::pair<unsigned int, Real>
TiberLinearSolver::solve(const ShellMatrix<Number>&,
    NumericVector<Number>&,
    NumericVector<Number>&,
    const double,
    const unsigned int)
{
  Messages::error("Solving with shell matrix is not implemented");
}
  


std::pair<unsigned int, Real>
TiberLinearSolver::solve(const ShellMatrix<Number>& shell_matrix,
    const SparseMatrix<Number>& precond_matrix,
    NumericVector<Number>&,
    NumericVector<Number>&,
    const double,
    const unsigned int)
{
  Messages::error("Solving with shell matrix is not implemented");
}
  


