// $Id$


#include "TiberNonlinearSystem.h"
#include "TiberPetscLinearSolver.h"
#include "InitFailedException.h"

#include "PetscDivergedError.h"


#include "linear_solver.h"
#include "equation_systems.h"


using namespace std;


TiberNonlinearSystem::TiberNonlinearSystem(EquationSystems& es,
    const string& name, const unsigned int number)
: Parent(es, name, number),
  _n_nonlin_iterations(0),
  _final_residual_norm(1e20),
  _last_step_size(1e20),
  _nonlin_step_tol(1e-3),
  _nonlin_rel_tol(1e-15),
  _nonlin_abs_tol(1e-15),
  _nonlin_max_it(20),
  _lin_tol(1e-6),
  _lin_max_it(500)
{
  _solver = new TiberPetscLinearSolver<double>();

  if (_solver == NULL)
    throw InitFailedException("Cannot create linear solver object.");

  // add a vector for the solution
  add_vector("sol");
}


TiberNonlinearSystem::~TiberNonlinearSystem(void)
{
  clear();

  delete _solver;
}




void
TiberNonlinearSystem::reinit(void)
{
  _solver->clear();

  Parent::reinit();
}



void
TiberNonlinearSystem::clear(void)
{
  _solver->clear();

  Parent::clear();
}



void
TiberNonlinearSystem::set_linear_solver_type(SolverType solver_type,
    PreconditionerType pc)
{
  _solver->set_solver_type(solver_type);
  _solver->set_preconditioner_type(pc);
}




void
TiberNonlinearSystem::solve(void)
{

  assert(_assemble != NULL);

  EquationSystems& es = get_equation_systems();

  NumericVector<Number>& u = get_vector("sol");
  NumericVector<Number>& du = *solution;
  AutoPtr<NumericVector<Number> > u_old_ptr = u.clone();
  NumericVector<Number>& u_old = *u_old_ptr;
  
  // the l_infty tolerance for the step size
  double eps = _nonlin_step_tol;
  
  // the tolerance for the residual
  double eps_res = _nonlin_abs_tol;

  double tol = _lin_tol;

  int max_ls_step = 10;

  bool cauchy = false;

  // the (final) residual norm
  double norm_res, norm_rhs = 0;

  // the norm of the search step
  double norm_du = 1e12, norm_du_old, norm_du_first = 1e12; 

  unsigned int i = 1;
  for ( ; i <= _nonlin_max_it; i++)
  {

    // prepare jacobian and residual
    _assemble(u, rhs, NULL);
    _assemble(u, NULL, matrix);
    
    // solve the linear system
    _solver->solve(*matrix, *solution, *rhs, tol, _lin_max_it);

    // the l2 norm of the current residual
    norm_rhs = rhs->l2_norm();
    norm_res = norm_rhs;

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

      cerr << "       ||r(x)|| = " << norm_rhs << ", ||r(x + " << 
        alpha << "*dx)|| = " << norm_res << endl;

      if ((norm_res < norm_rhs) || (norm_du < eps) || (norm_res < eps_res))
        break;
      else
      {
        if (ls_step == max_ls_step - 1)
        {
          //cauchy = true;

          //u.add(-1.0, du);
          //break;
        }

        // don't accept step
        u = u_old;
        alpha *= 0.5;

        // if the norm decreases sufficiently, we don't increase the counter
        if (norm_res < 0.5 * old_norm)
          ls_step--;
      }
    }

    if ((norm_du < eps) || (norm_res < eps_res))
      break;

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
      cerr << "        ||r(x + " << alpha << "*dx)|| = " << norm_res_old <<
        ", ||r(x + " << 0.5 * alpha << "*dx)|| = "  << norm_res << endl;
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


    cerr << "  it " << i << ", |du| = " << norm_du << ", |r| = " << norm_res << endl;

    tol *= tol;
      

    //if (norm_du < eps)
    if ((norm_du < eps) || (norm_res < eps_res))
      break;


    //ostringstream o;
    //o << "du_" << i << ".gmv";
    //vector<string> names;
    //vector<Number> sol;
    //es.build_solution_vector(sol);
    //es.build_variable_names(names);
    //GraceIO(get_mesh()).write_nodal_data(o.str(), sol, names);

  }

  _n_nonlin_iterations = i;
  _final_residual_norm = norm_res;
  _last_step_size = norm_du;

  cerr << "iterations: " << i << ", |du| = " << norm_du
    << ", |r| = " << norm_res << endl;

  if (i >= _nonlin_max_it)
    throw (PetscDivergedError(-3, i, norm_rhs));


  update();
}
