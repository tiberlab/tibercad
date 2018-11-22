// $Id$


#include "TiberNonlinTR.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#include "InitFailedException.h"

#include "SolveFailedException.h"


#include "libmesh/linear_solver.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

#include "Messages.h"

#include <cassert>



using namespace std;


TiberNonlinTR::TiberNonlinTR(libMesh::EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number)
{
}



TiberNonlinTR::~TiberNonlinTR(void)
{
}





void
TiberNonlinTR::do_solve(void)
{

  assert(_assemble != NULL);

  // the ghosted solution vector
  NumericVector<Number>& u = get_solution_vector();
  if (!u.closed()) u.close();

  // the linear system solution
  libMesh::UniquePtr<NumericVector<Number> > du_ptr = solution->clone();
  NumericVector<Number>& du = *du_ptr;

  libMesh::UniquePtr<NumericVector<Number> > u_old_ptr = solution->clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  NumericVector<Number>& u_tmp = *solution;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();

  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  double delta_max = get_max_abs_step(); // > 0
  double delta =  delta_max / 5.0; // = (0, delta)
  double eta = 0.1; // = [0, 1/4)

  // the (final) residual norm
  double norm_rhs = 1e56;

  // the norm of the search step
  double norm_du = 1e12;

  double tol = get_linear_solver()->get_linear_rtol();
  double tol_orig = tol;

  // for testing
  //UniquePtr<SparseMatrix<Number> > transpose;
  libMesh::UniquePtr<NumericVector<double> > gradf = rhs->clone();
  libMesh::UniquePtr<NumericVector<double> > tmpvec = solution->clone();

  // this is only when l2 norm is used
  //double sqrtn = std::sqrt(du.size());
  //delta_max *= sqrtn;
  //delta *= sqrtn;

  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    // prepare jacobian and residual
    _assemble(u, NULL, matrix, *this);
    _assemble(u, rhs, NULL, *this);

    double norm_rhs_now = rhs->l2_norm();

    get_linear_solver()->set_linear_rtol(tol);

    //
    // calculate cauchy point
    //
    matrix->get_transpose(*matrix); // J'
    matrix->vector_mult(*gradf, *rhs);
    double norm_pk_l2 = gradf->l2_norm();
    double norm_pk_infty = gradf->linfty_norm();
    //cerr << "|Z| = " << norm_rhs_now << "  " << " |g| = " << norm_grad_l2 << endl;

    //cerr << "|f| = " << norm_rhs_now << "  |g| = " << norm_grad_l2 << endl;

    // r'J(J'J)J'r
    matrix->get_transpose(*matrix); // J
    matrix->vector_mult(du, *gradf); // JJ'r
    matrix->get_transpose(*matrix);  // J'
    matrix->vector_mult(*tmpvec, du); //J'JJ'r
    double tmp = gradf->dot(*tmpvec); // r'JJ'JJ'r
    //double fac = tmp / (norm_grad_l2 * norm_grad_l2);
    //cerr << "g'(J'J)g = " << tmp << " fac = "<< fac << endl;
    tmp *= delta; // delta * r'JJ'JJ'r
    tmp = norm_pk_l2 * norm_pk_l2 * norm_pk_infty / tmp;
    double tau = std::min(1.0, tmp);
    gradf->scale(-tau * delta / norm_pk_infty);  // Cauchy point

    matrix->get_transpose(*matrix);  // J

    // length of the Cauchy step
    norm_du = gradf->linfty_norm();
    //norm_du = gradf->l2_norm();

    //cerr << "tau = " << tau << " |du| = " << norm_du << " delta = " << delta << endl;

    if (tau == 1.0)
      // we can just take the Cauchy step
      du = *gradf;
    else
    {
      //
      // calculate unconstrained minimizer (Newton step)
      //
      get_linear_solver()->solve(*matrix, du, *rhs);
      du.scale(-1.0);
      du.close();

      double t = 1.0;
      norm_du = du.linfty_norm();
      //double norm_newton = du.l2_norm();
      //cerr << "ratio = " << (norm_newton / norm_cauchy) << endl;
      //cerr << "  |du| = " << norm_du << " (delta = " << delta << ")\n";

      if (norm_du > delta)
      {
        double t0 = 0, t1 = 1;
        int i = 0, imax = 7;
        while (i < imax)
        {
          t = 0.5 * (t0 + t1);

          *tmpvec = du;
          tmpvec->scale(t);
          tmpvec->add(1 - t, *gradf);

          norm_du = tmpvec->linfty_norm();
          //norm_du = tmpvec->l2_norm();
          //cerr << "----- t = " << t << " |du| = " << norm_du << endl;

          if (norm_du > delta)
            t1 = t;
          else
            t0 = t;

          ++i;
        }

        // to be sure that we are inside the trust radius
        if (norm_du > delta)
        {
          //t1 = t;
          //t = 0.5 * (t0 + t1);
          t = t0;
          *tmpvec = du;
          tmpvec->scale(t);
          tmpvec->add(1 - t, *gradf);

          norm_du = tmpvec->linfty_norm();
        }

        du = *tmpvec;
      }
      //cerr << "calculated pc: t = " << t << " |du| = " << norm_du << endl;
    }


    // now we have the step p_k

    //
    // calculate rho_k = actual reduction / predicted reduction
    //

    u_tmp = u;
    u_old = u;
    u_tmp.add(du);
    u_tmp.localize(u, get_dof_map().get_send_list());


    //if ((norm_du < eps) || (norm_rhs < eps_res))
    if (norm_du < eps)
    {
      break;
    }
    else if (i == get_nonlinear_max_it())
    {
      throw (SolveFailedException("Trust-region method failed."));
    }


    _assemble(u, tmpvec.get(), NULL, *this);

    norm_rhs = tmpvec->l2_norm();

    matrix->vector_mult(*tmpvec, du);
    rhs->add(*tmpvec);
    double rhs_norm_pred = rhs->l2_norm();

    double r = norm_rhs_now * norm_rhs_now;
    double rho = (r - norm_rhs * norm_rhs) / (r - rhs_norm_pred * rhs_norm_pred);

    //cerr << "real = " << (r - norm_rhs * norm_rhs) <<
    //    " pred = " << (r - rhs_norm_pred * rhs_norm_pred) << "  rho = " << rho << endl;

    double delta_old = delta;

    if (rho < 0.25)
      delta = 0.25 * norm_du;
    else
    {
      if ((rho > 0.75) && Utils::almost_equal::compare(norm_du, delta, 1e-16))
        delta = std::min(2 * delta, delta_max);
      //else delta_(k+1) = delta_k
    }


    if (rho < eta)
      u = u_old;
    else
    {
      //u += du

      //norm_du = du.linfty_norm();
      {
        ostringstream os;
        os << "it " << i << ", |du| = " << norm_du
        << ", |r| = " << norm_rhs << " delta = " << delta_old;
        Messages::info(os.str());
      }

      draw_point(i, norm_du);

      tol *= tol;
    }


    // check for divergence
    //if ((norm_res > norm_rhs) || isnan(norm_res))
    //if (isnan(norm_res))
    //{
      //cout << endl;
    //  throw (SNESDivergedError(-4, i, norm_rhs));
    //}


    tol *= tol;

  }

  // reset the original tolerance
  get_linear_solver()->set_linear_rtol(tol_orig);

  _n_nonlin_iterations = i;
  _final_residual_norm = norm_rhs;
  _last_step_size = norm_du;

  ostringstream os;
  os << "iterations: " << i << ", |du| = " << norm_du
    << ", |r| = " << norm_rhs << Messages::endl;
  Messages::newline();
  Messages::info(os.str());


  update();
}
