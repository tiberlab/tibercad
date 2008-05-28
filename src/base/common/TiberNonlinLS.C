// $Id$


#include "TiberNonlinLS.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#include "InitFailedException.h"

#include "PetscDivergedError.h"
#include "SNESDivergedError.h"


#include "linear_solver.h"
#include "equation_systems.h"

#define DEBUG

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
  NumericVector<Number>& du = *solution;
  AutoPtr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  // the l_infty tolerance for the step size
  double eps = get_nonlinear_stol();
  
  // the tolerance for the residual
  double eps_res = get_nonlinear_atol();


  int max_ls_step = 5;

  bool cauchy = false;

  // the (final) residual norm
  double norm_res, norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12, norm_du_old; 

  // let the solver know the options
  get_linear_solver()->set_options(get_options());
  double tol = get_linear_solver()->get_linear_rtol();


  unsigned int i = 1;
  for ( ; i <= get_nonlinear_max_it(); i++)
  {

    // prepare jacobian and residual
    _assemble(u, rhs, NULL);
    _assemble(u, NULL, matrix);

    get_linear_solver()->set_linear_rtol(tol);
    
    // solve the linear system
    get_linear_solver()->solve(*matrix, *solution, *rhs);
#ifndef DEBUG
    cout << "." << flush;
#endif

    // the l2 norm of the current residual
    norm_rhs = rhs->l2_norm();
    norm_res = norm_rhs;
    
    if (norm_res < eps_res)
    {
#ifndef DEBUG
      cout << endl;
#endif
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
      u.add(-alpha, du);

      // evaluate the residual
      _assemble(u, rhs, NULL);

      old_norm = norm_res;
      norm_res = rhs->l2_norm();
      norm_du = du.linfty_norm();

      //cerr << "       ||r(x)|| = " << norm_rhs << ", ||r(x + " << 
      //  alpha << "*dx)|| = " << norm_res << endl;

      if ((norm_res < norm_rhs) || (norm_du < eps) || (norm_res < eps_res))
        break;
      else
      {
        if (ls_step == max_ls_step - 1)
        {
          //cauchy = true;

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
#ifndef DEBUG
      cout << endl;
#endif
      break;
    }
      

    // check for divergence
    //if ((norm_res > norm_rhs) || isnan(norm_res))
    if (isnan(norm_res))
    {
#ifndef DEBUG
      cout << endl;
#endif
      throw (SNESDivergedError(-4, i, norm_rhs));
    }

  

    if (cauchy)
    //if ((ls_step >= max_ls_step) && (norm_res >= norm_rhs))
    {
      cerr << "  trying Cauchy step:" << endl;
      cauchy = false;
      
      //u = u_old;
      //u.add(-1.0, du);
      //continue;
      
      //u.add(-0.5, du);
      u.add(-0.1 / rhs->l2_norm(), *rhs);

      // evaluate the residual
      _assemble(u, rhs, NULL);

      norm_res = rhs->l2_norm();
      norm_du = du.linfty_norm();
    }
    else
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

    //if (norm_du > _max_step_size)
    //{
    //  du.scale(_max_step_size / norm_du);
    //  norm_du = _max_step_size;
    //}

#ifdef DEBUG
    cout << "  it " << i << ", |du| = " << norm_du << ", |r| = " << norm_res << endl;
#endif

    draw_point(i, norm_res);

    tol *= tol;
      

    //if (norm_du < eps)
    if ((norm_du < eps) || (norm_res < eps_res))
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
  _final_residual_norm = norm_res;
  _last_step_size = norm_du;

  cout << "iterations: " << i << ", |du| = " << norm_du
    << ", |r| = " << norm_res << endl;

  
  update();
}
