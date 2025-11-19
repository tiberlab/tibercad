// $Id$

#include "tibercad/solver/TiberLineSearch.h"
#include "core/solver-interfaces/TiberLinearSolver.h"
#include "tibercad/base/InitFailedException.h"

#include "tibercad/base/TiberCad.h"

#include "mesh.h"

using namespace std;


TiberLineSearch::TiberLineSearch(libMesh::EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number),
  _solver(NULL),
  _nonlinear_rtol(1e-6),
  _nonlinear_atol(1e-50),
  _nonlinear_stol(1e-3),
  _nonlinear_max_it(25),
  _max_abs_step(1e3),
  _max_step(1e3),
  _divergence_tol(4.0)
{
  // add a vector for the solution
  //add_vector("sol", true, libMesh::GHOSTED);
}


TiberLineSearch::~TiberLineSearch(void)
{  
  if (_solver != NULL)
    _solver->clear();
  
  delete _solver;
}



void
TiberLineSearch::reinit(void)
{
  _solver->clear();

  // we could have changed the solver type
  user_initialization();

  // let the linear solver know the options
  ModelOptions::submodel_iterator it = get_options().submodels_begin("linear_solver");
  if (it != get_options().submodels_end("linear_solver"))
    get_linear_solver()->set_options(it->second);


  Parent::reinit();
}



void
TiberLineSearch::clear(void)
{
  if (_solver != NULL)
    _solver->clear();

  Parent::clear();
}




void
TiberLineSearch::user_initialization(void)
{
   
  if (_solver == NULL)
  {
    ModelOptions::submodel_iterator it = get_options().submodels_begin("linear_solver");
    if (it != get_options().submodels_end("linear_solver"))
      _solver = TiberLinearSolver::create(this->comm(), it->second);
    else
      _solver = TiberLinearSolver::create(this->comm(), ModelOptions());
  }
}


void
TiberLineSearch::parse_options(void)
{
  _nonlinear_rtol = get_options().get_option("relative_tolerance", _nonlinear_rtol);
  _nonlinear_atol = get_options().get_option("absolute_tolerance", _nonlinear_atol);
  _nonlinear_stol = get_options().get_option("step_tolerance", _nonlinear_stol);
  _nonlinear_max_it = get_options().get_option("max_iterations", _nonlinear_max_it);
  _divergence_tol = get_options().get_option("divergence_tolerance", _divergence_tol);

  // setup the max line search step
  double sqrt_nn = std::sqrt((double) get_mesh().n_nodes() * n_vars());
  _max_abs_step = get_options().get_option("max_step", 10.0);
  _max_step = _max_abs_step * sqrt_nn;
}
