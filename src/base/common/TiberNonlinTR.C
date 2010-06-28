// $Id$


#include "TiberNonlinTR.h"
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


TiberNonlinTR::TiberNonlinTR(EquationSystems& es,
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

  NumericVector<Number>& u = get_solution_vector();
  if (!u.closed()) u.close();
  NumericVector<Number>& du = *solution;
  AutoPtr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();
  
  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  //int max_ls_step = 5;

  double delta_max = 10; // > 0
  double delta = 1; // = (0, delta)
  double eta = 0.1; // = [0, 1/4)

  // the (final) residual norm
  double norm_res, norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12, norm_du_old; 

  // let the solver know the options
  get_linear_solver()->set_options(get_options());
  double tol = get_linear_solver()->get_linear_rtol();

  // for testing
  //AutoPtr<SparseMatrix<Number> > transpose;
  AutoPtr<NumericVector<double> > gradf = rhs->clone();
  AutoPtr<NumericVector<double> > tmpvec = solution->clone();

  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    // prepare jacobian and residual
    _assemble(u, NULL, matrix);
    _assemble(u, rhs, NULL);

    get_linear_solver()->set_linear_rtol(tol);
    
    //
    // calculate cauchy point
    //
    matrix->get_transpose(*matrix);
    matrix->vector_mult(*gradf, *rhs);
    double norm_grad = gradf->l2_norm();

    // r'J(J'J)J'r
    du = *gradf;
    matrix->get_transpose(*matrix); // J
    matrix->vector_mult(du, *gradf); // JJ'r
    matrix->get_transpose(*matrix);  // J'
    matrix->vector_mult(*tmpvec, du); //J'JJ'r
    double tmp = delta * gradf->dot(*tmpvec); // r'JJ'JJ'r
    tmp = norm_grad * norm_grad * norm_grad / tmp;
    double tau = std::min(1.0, tmp);
    gradf->scale(-tau * delta / norm_grad);  // Cauchy point

    matrix->get_transpose(*matrix);  // J

    norm_du = gradf->l2_norm();

    cerr << "tau = " << tau << " |du| = " << norm_du << endl;

    if (Utils::almost_equal::compare(norm_du, delta))
      du = *gradf;
    else
    {
      //
      // calculate unconstrained minimizer (Newton step)
      //
      get_linear_solver()->solve(*matrix, *solution, *rhs);
      du.scale(-1.0);

      double t = 1.0;
      norm_du = du.l2_norm();
      cerr << "  |du| = " << norm_du << " (delta = " << delta << ")\n";

      if (norm_du > delta)
      {
        double t0 = 0, t1 = 1;
        int i = 0, imax = 5;
        while (i < imax)
        {
          double t = 0.5 * (t0 + t1);

          *tmpvec = du;
          tmpvec->scale(t);
          tmpvec->add(1 - t, *gradf);

          norm_du = tmpvec->l2_norm();

          if (norm_du > delta)
            t1 = t;
          else
            t0 = t;

          ++i;
        }

        du = *tmpvec;
      }
      cerr << "calculated pc: t = " << t << " |du| = " << norm_du << endl;
    }


    // now we have the step p_k

    //
    // calculate rho_k = actual reduction / predicted reduction
    //
    double norm_rhs_now = rhs->l2_norm();

    u_old = u;
    u.add(du);
    _assemble(u, tmpvec.get(), NULL);

    norm_rhs = tmpvec->l2_norm();

    matrix->vector_mult(*tmpvec, du);
    rhs->add(*tmpvec);
    double rhs_norm_pred = rhs->l2_norm();

    double r = norm_rhs_now * norm_rhs_now;
    double rho = (r - norm_rhs * norm_rhs) / (r - rhs_norm_pred * rhs_norm_pred);

    cerr << "real = " << (r - norm_rhs * norm_rhs) <<
        " pred = " << (r - rhs_norm_pred * rhs_norm_pred) << "  rho = " << rho << endl;

    if (rho < 0.25)
      delta = 0.25 * norm_du;
    else
    {
      if ((rho > 0.75) && Utils::almost_equal::compare(norm_du, delta))
        delta = std::min(2 * delta, delta_max);
      //else delta_(k+1) = delta_k
    }
    
    if (rho < eta)
      u = u_old;
    // else u += du



    // check for divergence
    //if ((norm_res > norm_rhs) || isnan(norm_res))
    //if (isnan(norm_res))
    //{
      //cout << endl;
    //  throw (SNESDivergedError(-4, i, norm_rhs));
    //}

    /*
    {
      ostringstream os;
      os << "it " << i << ", |du| = " << norm_du
        << ", |r| = " << norm_res;
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
    */

  }

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
