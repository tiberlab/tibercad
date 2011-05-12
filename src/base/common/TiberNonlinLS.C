// $Id$


#include "TiberNonlinLS.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#include "InitFailedException.h"

#include "SolveFailedException.h"
#include "PetscDivergedError.h"
#include "SNESDivergedError.h"


#include "linear_solver.h"
#include "equation_systems.h"

#include "Messages.h"

#include <cassert>



using namespace std;


TiberNonlinLS::TiberNonlinLS(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number)
{
}



TiberNonlinLS::~TiberNonlinLS(void)
{
}





void
TiberNonlinLS::do_solve(void)
{

  assert(_assemble != NULL);

  NumericVector<Number>& u = get_solution_vector();
  if (!u.closed()) u.close();
  NumericVector<Number>& du = *solution;
  AutoPtr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();

  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  int max_ls_step = 5;

  // the (final) residual norm
  double norm_res, norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12, norm_du_old;

  double tol = get_linear_solver()->get_linear_rtol();
  double tol_orig = tol;


  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    // prepare jacobian and residual
    _assemble(u, NULL, matrix);
    _assemble(u, rhs, NULL);

    get_linear_solver()->set_linear_rtol(tol);

    // solve the linear system
    get_linear_solver()->solve(*matrix, *solution, *rhs);

    //cout << "." << flush;

    //matrix->vector_mult(*rhs, *solution);
    //double d = rhs->dot(*rhs);
    //cerr << "d = " << d << endl;

    // the l2 norm of the current residual
    norm_rhs = rhs->l2_norm();
    norm_res = norm_rhs;

    if (norm_res < eps_res)
    {
      //cout << endl;
      norm_du = 0.0;
      break;
    }

    u_old = u;
    norm_du_old = norm_du;

    // the relaxation factor
    double alpha = 1.0;

    /*
     * Note about the LS algorithm
     * ---------------------------
     *
     * Somehow we have to fix the number of line-search steps. We don't do
     * this, however, using a fixed maximum number of LS steps. If the
     * l2 norm of the residual between two LS steps decreases by more than
     * 50%, we don't increase the step counter.
     *
     */
    double old_norm;
    int ls_step = 0;
    for ( ; ls_step < max_ls_step; ls_step++)
    {
      // apply step and look at the new residual

      //norm_du = du.linfty_norm();
      //if (norm_du > get_max_abs_step())
      //{
      //  double fac = get_max_abs_step() / norm_du;
      //  du.scale(fac);
      //  norm_du *= fac;
      //}

      u.add(-alpha, du);

      // evaluate the residual
      _assemble(u, rhs, NULL);

      old_norm = norm_res;
      norm_res = rhs->l2_norm();
      norm_du = du.linfty_norm();

      // TODO this seems not to be a brilliant idea
      //if (norm_du > get_divergence_tol() * norm_du_old)
      //  throw (SolveFailedException("Line search diverged"));

      //cerr << "       ||r(x)|| = " << norm_rhs << ", ||r(x + " <<
      //  alpha << "*dx)|| = " << norm_res << endl;


      if ((norm_res < norm_rhs) || (norm_du < eps) || (norm_res < eps_res))
        break;
      else
      {
        if (ls_step == max_ls_step - 1)
        {
          //u.add(-1.0, du);
          break;
        }

        // don't accept step
        u = u_old;
        alpha *= 0.5;

        // if the norm decreases sufficiently, we don't increase the counter
        if (norm_res < 0.5 * old_norm)
          ls_step--;
      }
    }

    // check for convergence
    if ((norm_du < eps) || (norm_res < eps_res))
    {
      //cout << endl;
      break;
    }


    // check for divergence
    //if ((norm_res > norm_rhs) || isnan(norm_res))
    if (std::isnan(norm_res))
    {
      //cout << endl;
      throw (SNESDivergedError(-4, i, norm_rhs));
    }

    ///*
    {
      // check for one more smaller step
      u = u_old;
      u.add(-0.5 * alpha, du);

      // evaluate the residual
      _assemble(u, rhs, NULL);

      double norm_res_old = norm_res;
      norm_res = rhs->l2_norm();
      //cerr << "        ||r(x + " << alpha << "*dx)|| = " << norm_res_old <<
      //  ", ||r(x + " << 0.5 * alpha << "*dx)|| = "  << norm_res << endl;
      if (norm_res > norm_res_old)
      {
        // keep the former step
        u = u_old;
        u.add(-alpha, du);
        norm_res = norm_res_old;
      }
      else
        alpha *= 0.5;

      du.scale(alpha);
      norm_du *= alpha;
    }
    //*/

    //if (norm_du > _max_step_size)
    //{
    //  du.scale(_max_step_size / norm_du);
    //  norm_du = _max_step_size;
    //}

    {
      ostringstream os;
      os << "it " << i << ", |du| = " << norm_du
        << ", |r| = " << norm_res << "  alpha = " << alpha;
      Messages::info(os.str());
    }

    draw_point(i, norm_res);

    tol *= tol;


    //if (norm_du < eps)
    if ((norm_du < eps) || (norm_res < eps_res))
    {
      //cout << endl;
      break;
    }
    else if (i == get_nonlinear_max_it())
    {
      //cout << endl << flush;
      throw (PetscDivergedError(-3, i, norm_rhs));
    }

  }

  // reset the old linear tolerance
  get_linear_solver()->set_linear_rtol(tol_orig);

  _n_nonlin_iterations = i;
  _final_residual_norm = norm_res;
  _last_step_size = norm_du;

  ostringstream os;
  os << "iterations: " << i << ", |du| = " << norm_du
    << ", |r| = " << norm_res << Messages::endl;
  Messages::newline();
  Messages::info(os.str());


  update();
}
