// $Id$


#include "TiberNonlinBR.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#include "InitFailedException.h"

#include "PetscDivergedError.h"
#include "SNESDivergedError.h"


#include "linear_solver.h"
#include "equation_systems.h"

#define DEBUG

using namespace std;


TiberNonlinBR::TiberNonlinBR(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number),
  _solver(NULL)
{
  // add a vector for the solution
  add_vector("sol");
}



TiberNonlinBR::~TiberNonlinBR(void)
{
  clear();

  delete _solver;
}




void
TiberNonlinBR::reinit(void)
{
  _solver->clear();

  _solver->set_solver_type(_solver_type);
  _solver->set_preconditioner_type(_preconditioner_type);

  // we could have changed the solver type
  user_initialization();

  Parent::reinit();
}



void
TiberNonlinBR::clear(void)
{
  if (_solver != NULL)
  {
    _solver->clear();

    _solver->set_solver_type(_solver_type);
    _solver->set_preconditioner_type(_preconditioner_type);
  }

  Parent::clear();
}



void
TiberNonlinBR::user_initialization(void)
{
  if (_solver == NULL)
  {
    _solver = TiberLinearSolver::create(_linear_solver);

    if (_solver == NULL)
    {
      std::cerr << "Linear solver " << _linear_solver <<
        " is not available." << std::endl;
      throw InitFailedException("Cannot create linear solver object.");
    }
  }

  _solver->set_solver_type(_solver_type);
  _solver->set_preconditioner_type(_preconditioner_type);
}


void
TiberNonlinBR::solve(void)
{

  assert(_assemble != NULL);

  NumericVector<Number>& u = get_vector("sol");
  NumericVector<Number>& du = *solution;
  AutoPtr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;

  // the l_infty tolerance for the step size
  double eps = _nonlin_step_tol;
  
  // the tolerance for the residual
  double eps_res = _nonlin_abs_tol;

  double tol = _lin_tol;


  // the (final) residual norm
  double norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12; 


  double d = 0.8;
  double K = 0;

  unsigned int i = 1;
  for ( ; i <= _nonlin_max_it; i++)
  {

    // prepare jacobian and residual
    _assemble(u, rhs, NULL);
    _assemble(u, NULL, matrix);

    _solver->set_ksp_options(tol, _lin_abs_tol, _lin_max_it);
    
    // solve the linear system
    _solver->solve(*matrix, *solution, *rhs, tol, _lin_max_it);
#ifndef DEBUG
    cout << "." << flush;
#endif

    // the l2 norm of the current residual
    norm_rhs = rhs->l2_norm();

    u_old = u;

    while (1)
    {
      double tk = 1.0 / (1.0 + K * norm_rhs);
    
      u.add(-tk, du);

      _assemble(u, rhs, NULL);
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

    
    norm_du = du.linfty_norm();


    // check for divergence
    if (isnan(norm_rhs))
    {
#ifndef DEBUG
      cout << endl;
#endif
      throw (SNESDivergedError(-4, i, norm_rhs));
    }


#ifdef DEBUG
    cout << "  it " << i << ", |du| = " << norm_du << ", |r| = " << norm_rhs << endl;
#endif

    // check for convergence
    if ((norm_du < eps) || (norm_rhs < eps_res))
    {
#ifndef DEBUG
      cout << endl;
#endif
      break;
    }
    else if (i == _nonlin_max_it)
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
