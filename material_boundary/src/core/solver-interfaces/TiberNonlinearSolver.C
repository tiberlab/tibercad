// $Id$

#include "TiberNonlinearSolver.h"
#include "TiberPetscNonlinearSolver.h"
#include "ModelOptions.h"
#include "XMonitor.h"

// default values
namespace
{
  const double default_nonlinear_rtol = 1e-9;
  const double default_nonlinear_atol = 1e-50;
  const double default_nonlinear_stol = 1e-3;
  const int default_nonlinear_max_it = 25;
  const double default_linear_rtol = 1e-6;
  const double default_linear_atol = 1e-50;
  const int default_linear_max_it = 500;
}


TiberNonlinearSolver::TiberNonlinearSolver(void)
  : _nonlinear_rtol(default_nonlinear_rtol),
    _nonlinear_atol(default_nonlinear_atol),
    _nonlinear_stol(default_nonlinear_stol),
    _nonlinear_max_it(default_nonlinear_max_it),
    _linear_rtol(default_linear_rtol),
    _linear_atol(default_linear_atol),
    _linear_max_it(default_linear_max_it),
    _xmonitor(NULL)
{
}



TiberNonlinearSolver*
TiberNonlinearSolver::create(const ModelOptions& options)
{
  return new TiberPetscNonlinearSolver();
}


void
TiberNonlinearSolver::set_options(const ModelOptions& options)
{
  _nonlinear_rtol = options.get_option("nonlin_rel_tol", default_nonlinear_rtol);
  _nonlinear_atol = options.get_option("nonlin_abs_tol", default_nonlinear_atol);
  _nonlinear_max_it = options.get_option("nonlin_max_it", default_nonlinear_max_it);
  _nonlinear_stol = options.get_option("nonlin_step_tol", default_nonlinear_stol);

  _linear_rtol = options.get_option("lin_rel_tol", default_linear_rtol);
  _linear_atol = options.get_option("lin_abs_tol", default_linear_atol);
  _linear_max_it = options.get_option("lin_max_it", default_linear_max_it);

  parse_options(options);
}



void
TiberNonlinearSolver::draw_point(double iteration, double error, bool logarithm)
{
  if (_xmonitor != NULL)
  {
    if (logarithm)
      error = log10(error);

    _xmonitor->draw_point(iteration, error);
  }
}
