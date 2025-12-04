/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TiberNonlinearSolver.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "TiberNonlinearSolver.h"
#include "petsc/TiberPetscNonlinearSolver.h"
#include "tibercad/base/ModelOptions.h"
#include "XMonitor.h"

// default values
namespace
{
  const double default_nonlinear_rtol = 1e-9;
  const double default_nonlinear_atol = 1e-50;
  const double default_nonlinear_stol = 1e-3;
  const int default_nonlinear_max_it = 25;
}


TiberNonlinearSolver::TiberNonlinearSolver(sys_type& s)
  : libMesh::NonlinearSolver<double>(s),
    _nonlinear_rtol(default_nonlinear_rtol),
    _nonlinear_atol(default_nonlinear_atol),
    _nonlinear_stol(default_nonlinear_stol),
    _nonlinear_max_it(default_nonlinear_max_it),
    _xmonitor(NULL)
{
}


/*
TiberNonlinearSolver*
TiberNonlinearSolver::create(const ModelOptions& options)
{
  return new TiberPetscNonlinearSolver();
}
*/


void
TiberNonlinearSolver::set_options(const ModelOptions& options)
{
  _nonlinear_rtol = options.get_option("relative_tolerance", default_nonlinear_rtol);
  _nonlinear_atol = options.get_option("absolute_tolerance", default_nonlinear_atol);
  _nonlinear_max_it = options.get_option("max_iterations", default_nonlinear_max_it);
  _nonlinear_stol = options.get_option("step_tolerance", default_nonlinear_stol);

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
