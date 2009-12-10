// $Id$

#include "RelaxationMethod.h"
#include "Control.h"
#include "Messages.h"


using namespace std;



void
RelaxationMethod::parse_options(void)
{
  SelfconsistentSolver::parse_options();

  ModelOptions& opts = get_options();

  _relax = opts.get_option("relaxation_factor", _relax);
}




void
RelaxationMethod::do_solve(void)
{
  parse_options();

  double relax = get_option("special_relaxation", 1.0);

  //initialize();

  open_xmonitor();

  AutoPtr<NumericVector<double> > x_old = NumericVector<double>::build();
  x_old->init(get_solution_vector());
  x_old->close();

  for (unsigned int it = 0; it < get_maximum_iterations(); it++)
  {
    *x_old = get_solution_vector();
    x_old->close();


    double x_old_norm = x_old->l2_norm();
    {
      //ostringstream file;
      //file << "X_" << it << ".m";
      //get_solution_vector().print_matlab(file.str());
    }

    get_last_simulation()->_relaxation = relax;
    relax = std::sqrt(relax);
    solve_simulations();

    get_solution_vector().close();
    get_solution_vector() -= *x_old;
    double norm = get_solution_vector().linfty_norm();
    double rel_err =  get_solution_vector().l2_norm() / x_old_norm;
    {
      //ostringstream file;
      //file << "dX_" << it << ".m";
      //get_solution_vector().print_matlab(file.str());
    }


    bool converged = true;

    if (get_monitor())
    {
      Messages m;
      ostringstream os;
      os << get_name() << " (Relaxation): iteration "
        << it << Messages::endl
        << "  correction (max norm):  " << norm << Messages::endl
        << "  relative error (l2)  :  " << rel_err;
      m.info(os.str());
    }

    draw_point(it, rel_err);

    // check for the difference between old and new solutions
    if ((norm > get_absolute_tolerance()) &&
        (rel_err > get_relative_tolerance()))
      converged = false;

    get_solution_vector().scale(_relax);
    get_solution_vector() += *x_old;

    if (converged)
      break;

  }

  close_xmonitor();
}



