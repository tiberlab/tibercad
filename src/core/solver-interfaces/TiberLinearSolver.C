// $Id$

#include "tiber_config.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"

#include "TiberCad.h"

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




TiberLinearSolver::TiberLinearSolver(const libMesh::Parallel::Communicator &comm_in, const ModelOptions& options)
  : TiberModelObject(options), libMesh::LinearSolver<Number>(comm_in), 
    _linear_rtol(default_linear_rtol),
    _linear_atol(default_linear_atol),
    _linear_max_it(default_linear_max_it)
{

}




TiberLinearSolver*
TiberLinearSolver::create(const libMesh::Parallel::Communicator &comm_in, const ModelOptions& options)
{
  TiberLinearSolver* solver = NULL;

  std::string type(options.get_name());
  if (type.empty())
    type = "petsc";
  type = options.get_option("type", type);

  if (type == "petsc")
    solver = new TiberPetscLinearSolver(comm_in, options); 
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

  solver->parse_options();

  return solver;
}




void
TiberLinearSolver::parse_options(void)
{
  _linear_rtol = get_option("relative_tolerance", default_linear_rtol);
  _linear_atol = get_option("absolute_tolerance", default_linear_atol);
  _linear_max_it = get_option("max_iterations", default_linear_max_it);

  do_parse_options();

  // dummy read
  get_option("type", "");
  get_option("simulation", "");

  get_options().check_unused();
}



std::string
TiberLinearSolver::get_simulation_name(void) const
{
  return get_option("simulation", "");
}


std::pair<unsigned int, Real>
TiberLinearSolver::solve(const libMesh::ShellMatrix<Number>&,
    NumericVector<Number>&,
    NumericVector<Number>&,
    const double,
    const unsigned int)
{
  Messages::error("Solving with shell matrix is not implemented");
  return std::make_pair(0, 0.0);
}
  


std::pair<unsigned int, Real>
TiberLinearSolver::solve(const libMesh::ShellMatrix<Number>&,
    const SparseMatrix<Number>&,
    NumericVector<Number>&,
    NumericVector<Number>&,
    const double,
    const unsigned int)
{
  Messages::error("Solving with shell matrix is not implemented");
  return std::make_pair(0, 0.0);
}



