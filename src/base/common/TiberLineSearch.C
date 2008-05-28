// $Id$

#include "TiberLineSearch.h"
#include "TiberLinearSolver.h"
#include "InitFailedException.h"

#include "mesh.h"

using namespace std;


TiberLineSearch::TiberLineSearch(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number),
  _solver(NULL),
  _nonlinear_rtol(1e-6),
  _nonlinear_atol(1e-50),
  _nonlinear_stol(1e-3),
  _nonlinear_max_it(25),
  _max_step(1e3)
{
  // add a vector for the solution
  add_vector("sol");
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
    _solver = TiberLinearSolver::create(get_options());
}


void
TiberLineSearch::parse_options(void)
{
  _nonlinear_rtol = get_options().get_option("nonlin_rel_tol", _nonlinear_rtol);
  _nonlinear_atol = get_options().get_option("nonlin_abs_tol", _nonlinear_atol);
  _nonlinear_stol = get_options().get_option("nonlin_step_tol", _nonlinear_stol);
  _nonlinear_max_it = get_options().get_option("nonlin_max_it", _nonlinear_max_it);

  // setup the max line search step
  double sqrt_nn = std::sqrt((double) get_mesh().n_nodes() * n_vars());
  double ls_max_step = get_options().get_option("ls_max_step", 10.0);
  _max_step = ls_max_step * sqrt_nn;
}
