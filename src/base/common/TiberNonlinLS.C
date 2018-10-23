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
#include "dof_map.h"

#include "Messages.h"

#include <cassert>



using namespace std;


TiberNonlinLS::TiberNonlinLS(libMesh::EquationSystems& es,
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

  // the ghosted solution vector
  NumericVector<Number>& u = get_solution_vector();
  if (!u.closed()) u.close();

  // the linear system solution
  libMesh::UniquePtr<NumericVector<Number> > du_ptr = solution->clone();
  NumericVector<Number>& du = *du_ptr;

  // the old local solution values
  libMesh::UniquePtr<NumericVector<Number> > u_old_ptr = solution->clone();
  NumericVector<Number>& u_old = *u_old_ptr;


  //libMesh::UniquePtr<NumericVector<Number> > u_tmp_ptr = u_old.clone();
  NumericVector<Number>& u_tmp = *solution;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();

  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  int max_ls_step = 10;

  // the (final) residual norm
  double norm_res, norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12, norm_du_old;

  double tol = get_linear_solver()->get_linear_rtol();
  double tol_orig = tol;

  //UniquePtr<SparseMatrix<double> > transpose = SparseMatrix<double>::build();

  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    u_old = u_tmp;

    // prepare jacobian and residual
    _assemble(u, NULL, matrix, *this);
    _assemble(u, rhs, NULL, *this);

    get_linear_solver()->set_linear_rtol(tol);

    du.zero();

    try {
      // solve the linear system
      if (this->have_matrix("Preconditioner"))
        get_linear_solver()->solve(*matrix, this->get_matrix("Preconditioner"),
            du, *rhs);
      else
        get_linear_solver()->solve(*matrix, du, *rhs);
    }
    catch (...)
    {
      // reset the old linear tolerance
      get_linear_solver()->set_linear_rtol(tol_orig);

      throw;
    }


    du.close();

    // the l2 norm of the current residual
    //norm_rhs = rhs->l2_norm();
    norm_rhs = TiberEqSystem::calculate_norm(rhs, l2_NORM);
    norm_res = norm_rhs;


    //matrix->vector_mult(*rhs, *solution);
    //matrix->get_transpose(*transpose);
    //transpose->vector_mult(*tmp_vec, *rhs);
    //double costheta = norm_rhs * norm_rhs / (solution->l2_norm() * tmp_vec->l2_norm());
    //cerr << "cos(theta) = " << costheta << endl;


    if (norm_res < eps_res)
    {
      //cout << endl;
      norm_du = 0.0;
      break;
    }

    norm_du_old = norm_du;

    // the relaxation factor
    double alpha = 1.0;
    double min_alpha = 1.0;
    //max_ls_step = 2;

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

      /* This could be useful in some cases, but we have to assure that
       * it does not limit the step size in other cases.
      norm_du = du.linfty_norm();
      if (norm_du > get_max_abs_step())
      {
        double fac = get_max_abs_step() / norm_du;
        du.scale(fac);
        norm_du *= fac;
      }
      */

      u_tmp.add(-alpha, du);
      u_tmp.localize(u, get_dof_map().get_send_list());

      // evaluate the residual
      _assemble(u, rhs, NULL, *this);

      old_norm = norm_res;
      //norm_res = rhs->l2_norm();
      norm_res = TiberEqSystem::calculate_norm(rhs, l2_NORM);
      norm_du = TiberEqSystem::calculate_norm(&du, MAX_NORM);
      //norm_du = du.linfty_norm();
      //norm_du = du.l2_norm();
      //double norm_u = u.l2_norm();
      //eps = get_nonlinear_stol() * norm_u + 1e-12;


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

        if ((ls_step == 1) && (old_norm < 100 * norm_rhs))
        {
          double b = 2*(norm_rhs + old_norm) - 4*norm_res;
          double a = old_norm - norm_rhs - b;

          min_alpha = -a / (2 * b);
          //cerr << "         min: alpha = " << min_alpha << endl;
          if ((min_alpha <= 1.0) && (min_alpha > 0))
            alpha = min_alpha;
        }
        else
          alpha *= 0.5;

        // don't accept step
        u_tmp = u_old;

        // if the norm decreases sufficiently, we don't increase the counter
        if (norm_res < 0.9 * old_norm)
          ls_step--;
      }
    }

    /*
    if (norm_res > old_norm)
    {
      alpha = 1;
      norm_res = old_norm;

      u_tmp = u_old;
      u_tmp.add(-alpha, du);
      u_tmp.localize(u, get_dof_map().get_send_list());
    }
    */

    // check for convergence
    if ((norm_du < eps) || (norm_res < eps_res))
    {
      //cout << endl;
      break;
    }


    // check for divergence
    if ((norm_res > norm_rhs) || std::isnan(norm_res))
    //if (std::isnan(norm_res))
    {
      // reset the old linear tolerance
      get_linear_solver()->set_linear_rtol(tol_orig);

      //cout << endl;
      throw (SNESDivergedError(-4, i, norm_rhs));
    }


    //if ((min_alpha <= 1.0) && (min_alpha > 0))
    if (10 * norm_res > norm_rhs)
    {
      // check for one more smaller step
      u_tmp = u_old;
      u_tmp.add(-0.5 * alpha, du);
      //u_tmp.add(-min_alpha, du);
      u_tmp.localize(u, get_dof_map().get_send_list());

      // evaluate the residual
      _assemble(u, rhs, NULL, *this);

      double norm_res_old = norm_res;
      //norm_res = rhs->l2_norm();
      norm_res = TiberEqSystem::calculate_norm(rhs, l2_NORM);
      //cerr << "        ||r(x + " << alpha << "*dx)|| = " << norm_res_old <<
      //  ", ||r(x + " << 0.5 * alpha << "*dx)|| = "  << norm_res << endl;
      if (norm_res > norm_res_old)
      {
        // keep the former step
        u_tmp = u_old;
        u_tmp.add(-alpha, du);
        u_tmp.localize(u, get_dof_map().get_send_list());
        norm_res = norm_res_old;
      }
      else
        alpha *= 0.5;
        //alpha = min_alpha;

      //du.scale(alpha);
      //norm_du *= alpha;
    }



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
    tol = (tol < 1e-12) ? 1e-12 : tol;


    //if (norm_du < eps)
    if ((norm_du < eps) || (norm_res < eps_res))
    {
      //cout << endl;
      break;
    }
    else if (i == get_nonlinear_max_it())
    {
      // reset the old linear tolerance
      get_linear_solver()->set_linear_rtol(tol_orig);

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
