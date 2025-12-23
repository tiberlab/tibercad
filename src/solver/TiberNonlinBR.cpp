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
 * \file TiberNonlinBR.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */



#include "TiberNonlinBR.h"
#include "solver/TiberLinearSolver.h"
#include "solver/petsc/TiberPetscLinearSolver.h"
#include "tibercad/base/InitFailedException.h"

#include "solver/petsc/PetscDivergedError.h"
#include "solver/petsc/SNESDivergedError.h"


#include "linear_solver.h"
#include "equation_systems.h"

#include <cassert>

#define DEBUG

using namespace std;


TiberNonlinBR::TiberNonlinBR(libMesh::EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number)
{
}



TiberNonlinBR::~TiberNonlinBR(void)
{
}



void
TiberNonlinBR::do_solve(void)
{

  assert(_assemble != NULL);

  //NumericVector<Number>& u = get_vector("sol");
  NumericVector<Number>& u = get_solution_vector();
  NumericVector<Number>& du = *solution;
  unique_ptr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();
  
  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  // the (final) residual norm
  double norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12; 


  double d = 0.8;
  double K = 0;

  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    // prepare jacobian and residual
    _assemble(u, rhs, NULL, *this);
    _assemble(u, NULL, matrix, *this);

    // solve the linear system
    get_linear_solver()->solve(*matrix, *solution, *rhs);

#ifndef DEBUG
    cout << "." << flush;
#endif

    // the l2 norm of the current residual
    norm_rhs = rhs->l2_norm();

    u_old = u;
    norm_du = du.linfty_norm();

    /*
    if (norm_du > get_max_abs_step())
    {
      double fac = get_max_abs_step() / norm_du;
      du.scale(fac);
      norm_du *= fac;
    }
    */

    /*
    double norm_du2 = du.l2_norm();
    if (norm_du2 > get_max_step())
    {
      du.scale(get_max_step() / norm_du2);
      norm_du = du.linfty_norm();
    }
    */

    while (1)
    {
      double tk = 1.0 / (1.0 + K * norm_rhs);
    
      u.add(-tk, du);

      _assemble(u, rhs, NULL, *this);
      double norm_rhs_new = rhs->l2_norm();

      // we check for convergence here to not get stuck in the inner loop
      if ((norm_du < eps) || (norm_rhs < eps_res))
        break;

      double e = (1.0 - norm_rhs_new / norm_rhs) / tk;
      if (e < d)
      {
        if (K == 0) K = 1;
        else K = 10 * K;

        u = u_old;
      }
      else
      {
        K = K / 10.0;

        // get out of the inner loop
        break;
      }
    }


    // check for divergence
    if (std::isnan(norm_rhs))
    {
#ifndef DEBUG
      cout << endl;
#endif
      throw (SNESDivergedError(-4, i, norm_rhs));
    }


#ifdef DEBUG
    cout << "  it " << i << ", |du| = " << norm_du << ", |r| = " << norm_rhs << endl;
#endif

    draw_point(i, norm_rhs);

    // check for convergence
    if ((norm_du < eps) || (norm_rhs < eps_res))
    {
#ifndef DEBUG
      cout << endl;
#endif
      break;
    }
    else if (i == get_nonlinear_max_it())
    {
#ifndef DEBUG
      cout << endl << flush;
#endif
      throw (PetscDivergedError(-3, i, norm_rhs));
    }

  }

  _n_nonlin_iterations = i;
  _final_residual_norm = norm_rhs;
  _last_step_size = norm_du;

  cout << "iterations: " << i << ", |du| = " << norm_du
    << ", |r| = " << norm_rhs << endl;

  
  update();
}
