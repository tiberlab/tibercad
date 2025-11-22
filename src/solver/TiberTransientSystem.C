#include "tibercad/solver/TiberTransientSystem.h"
#include "solver/TiberLinearSolver.h"
#include "tibercad/base/SolveFailedException.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/io/Messages.h"

#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"


template <class Base>
TiberTransientSystem<Base>::TiberTransientSystem(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number)
  :
    parent_type(es, name, number)
{
}


template <class Base>
TiberTransientSystem<Base>::~TiberTransientSystem(void) = default;


template <class Base>
TiberTransientSystem<Base>*
TiberTransientSystem<Base>::create(libMesh::EquationSystems& es,
        const std::string& sysname, const ModelOptions& options)
{
  TiberTransientSystem<Base>* sys = nullptr;

  sys = &(es.add_system<TiberTransientSystem<Base>>(sysname));

  if (sys == nullptr)
    throw InitFailedException("Could not create linear system " + sysname);

  sys->set_options(options);

  return sys;
}


template <class Base>
void
TiberTransientSystem<Base>::set_target_time(double time)
{
  _target_time = time;
}


template <class Base>
void
TiberTransientSystem<Base>::set_dirichlet_dofs(const std::set<libMesh::dof_id_type>& dirichlet_dofs)
{
  _dirichlet_dofs = dirichlet_dofs;
}


template <class Base>
void
TiberTransientSystem<Base>::solve(void)
{

  // if the time has not been updated, we can return back
  if (_target_time <= this->time)
    return;

  // update the old solution vector to the current one
  *this->old_local_solution = this->get_local_solution_vector();

  std::string method = "backward_euler";
  method = get_options().get_option("method", method);

  if (method == "backward_euler")
    backward_euler();
  else if (method == "forward_euler")
    forward_euler();
  else if (method == "trapezoidal")
    trapezoidal();
  else if (method == "adaptive")
    adaptive();
  else
  {
    Messages::error("Unknown transient solver: " + method);
    throw std::runtime_error("Unknown transient solver: " + method);
  }

}


template <class Base>
void
TiberTransientSystem<Base>::forward_euler(void)
{

  // the full time step to achieve
  double delta_t = _target_time - this->time;

  // set the systems time to the target, for assembly
  this->time = _target_time;

  NumericVector<libMesh::Number>& A = this->get_vector("t_weight");

  // I think this should be always true, here
  if (this->assemble_before_solve)
  {
    A.zero();
    this->assemble();
  }

  // Get a reference to the EquationSystems
  const libMesh::EquationSystems& es =
    this->get_equation_systems();


  // we calculate
  // 1/delta_t*A*c_new = 1/delta_t*A*c_old - K*c_old + f
  
  //*this->solution = *this->old_local_solution;

  // now we have basically a residual
  // we now set it explicitly to 0 for Dirichlet DoFs
  // We can figure out Dirichlet DoFs by looking at A, which
  // is assumed to be set to 0.0 at the corresponding indices
  for (auto i = A.first_local_index(); i < A.last_local_index(); ++i)
  {
    if (A.el(i) == 0.0)
    {
      A.set(i, delta_t);
      double val = this->rhs->el(i) / (*this->matrix)(i,i);
      this->solution->set(i, val);
    }
  }
  A.close();
  this->solution->close();

  // f -> f - K*c_old
  this->rhs->scale(-1.0);
  this->rhs->add_vector(*this->solution, *this->matrix);

  *this->rhs /= A;
  this->rhs->scale(-delta_t);
  
  this->solution->add(*this->rhs);

  // Update the system after the solve
  this->update();

}


template <class Base>
void
TiberTransientSystem<Base>::backward_euler(void)
{

  // the full time step to achieve
  double delta_t = _target_time - this->time;

  // set the systems time to the target, for assembly
  //double old_time = this->time;
  this->time = _target_time;

  NumericVector<libMesh::Number>& A = this->get_vector("t_weight");

  // I think this should be always true, here
  if (this->assemble_before_solve)
  {
    A.zero();
    this->assemble();
  }

  // Get a reference to the EquationSystems
  const libMesh::EquationSystems& es =
    this->get_equation_systems();

  TiberLinearSolver* lin_solver =
    static_cast<TiberLinearSolver*>(this->linear_solver.get());

  if (get_options().has_submodel("linear_solver"))
  {
    auto lin_opts = get_options().submodels_begin("linear_solver");
    lin_solver->set_options(lin_opts->second);
  }

  double lin_rel_tol = lin_solver->get_linear_rtol();
  int lin_max_it = lin_solver->get_linear_max_it();

  // setup the final system matrix
  // 1/delta_t*A + K
  // A is the diagonal, so we have to add element wise
  //
  // setup the final rhs
  // 1/delta_t*A*u_old + f
  //
  for (auto i = this->matrix->row_start(); i < this->matrix->row_stop(); ++i)
  {
    double mul = 1.0/delta_t * A.el(i);
    this->matrix->add(i, i, mul);
    this->rhs->add(i, mul * this->old_local_solution->el(i));
  }

  this->matrix->close();
  this->rhs->close();
  //this->matrix->print_matlab("K_solver.m");
  //this->rhs->print_matlab("F_solver.m");



  // Solve the linear system
  const std::pair<unsigned int, Real> rval;
  if (this->have_matrix("Preconditioner"))
    lin_solver->solve(*this->matrix, this->get_matrix("Preconditioner"),
        *this->solution, *this->rhs, lin_rel_tol, lin_max_it);
  else
    lin_solver->solve(*this->matrix, *this->solution, *this->rhs, lin_rel_tol, lin_max_it);


  // Store the number of linear iterations required to
  // solve and the final residual.
  this->_n_linear_iterations   = rval.first;
  this->_final_linear_residual = rval.second;

  // Update the system after the solve
  this->update();

  //this->time = _target_time;

  lin_solver->clear();
}


template <class Base>
void
TiberTransientSystem<Base>::adaptive(void)
{

  // the full time step to achieve
  double delta_t = _target_time - this->time;

  // the (last) relative error
  double rel_err = 0;

  // this is used to limit growing of step size
  // after a failed trial
  // the very first time we will not increase step size
  bool rejected = true;

  // reduce to the currently know working time step
  delta_t = std::min(delta_t, _time_step);

  NumericVector<libMesh::Number> &A = this->get_vector("t_weight");

  while (this->time < _target_time)
  {
    // adjust time, but if we overshoot the target, we have to 
    // adjust next_time_step, too
    if ((this->time + delta_t) > _target_time)
    {
      delta_t = _target_time - this->time;
    }
    this->time += delta_t;


    if (this->assemble_before_solve)
    {
      A.zero();
      this->assemble();
    }

    // We have now all matrices/vectors at the new time step.
    // We use this as a temporary simplification. In principle,
    // for time-dependent parameters, or in nonlinear case,
    // they should refer to the two different times.
   

    // Get a reference to the EquationSystems
    const libMesh::EquationSystems &es =
        this->get_equation_systems();

    TiberLinearSolver *lin_solver =
        static_cast<TiberLinearSolver *>(this->linear_solver.get());

    if (get_options().has_submodel("linear_solver"))
    {
      auto lin_opts = get_options().submodels_begin("linear_solver");
      lin_solver->set_options(lin_opts->second);
    }

    double lin_rel_tol = lin_solver->get_linear_rtol();
    int lin_max_it = lin_solver->get_linear_max_it();

    // set the boundary values
    for (auto i = A.first_local_index(); i < A.last_local_index(); ++i)
    {
      if (A.el(i) == 0.0)
      {
        double val = this->rhs->el(i) / (*this->matrix)(i, i);
        this->solution->set(i, val);
      }
    }
    this->solution->close();

    //
    // first, use trapezoidal method on the full step
    //

    auto tmp_vec = this->rhs->clone();

    // setup the final rhs for the trapezoidal method
    // (2/delta_t*A - K)*c_old + 2*f
    //
    this->solution->scale(-1.0);
    tmp_vec->scale(2.0);
    tmp_vec->add_vector(*this->solution, *this->matrix);
    (*this->solution) *= A;
    tmp_vec->add(-2.0 / delta_t, *this->solution);
    tmp_vec->close();

    // setup the final system matrix
    // 2/delta_t*A + K
    // A is the diagonal, so we have to add element wise
    //
    for (auto i = this->matrix->row_start(); i < this->matrix->row_stop(); ++i)
    {
      double val = 2.0 / delta_t * A.el(i);
      this->matrix->add(i, i, val);
    }
    this->matrix->close();

    // Solve the linear system, assume we have no preconditioner in this context
    std::pair<unsigned int, Real> rval;
    rval = lin_solver->solve(*this->matrix, *this->solution, *tmp_vec, lin_rel_tol, lin_max_it);

    //
    // now the backward Euler
    //

    // temporary solution vector used for the backward Euler
    auto backw_sol = this->solution->clone();

    for (auto i = this->matrix->row_start(); i < this->matrix->row_stop(); ++i)
    {
      double mul = 1.0 / delta_t * A.el(i);
      this->matrix->add(i, i, -mul);
      this->rhs->add(i, mul * this->old_local_solution->el(i));
    }
    this->matrix->close();
    this->rhs->close();

    rval = lin_solver->solve(*this->matrix, *backw_sol, *this->rhs, lin_rel_tol, lin_max_it);


    // now estimate the error
    backw_sol->add(-1.0, *this->solution);

    double err_norm = backw_sol->l2_norm();
    double l2_norm = this->solution->l2_norm() + 1e-12;
    rel_err = err_norm / l2_norm;
    
   

    double max_rel_err = get_options().get_option("relative_error", 1e-6);

    // limit change
    double max_inc = get_options().get_option("max_increase_factor", 5);
    double max_dec = get_options().get_option("min_decrease_factor", 0.2);

    // see Hairer/Nørsett/Wanner: "Solving ODE I" Chapter 4
    // maybe a different value could be tried also
    if (rejected)
      max_inc = 1.0;

    if (rel_err > max_rel_err)
    {
      this->time -= delta_t;
      *this->solution = *this->old_local_solution;
      rejected = true;
    }
    else
    {
      *this->old_local_solution = this->get_local_solution_vector();
      rejected = false;
    }

    double factor = std::sqrt(max_rel_err / rel_err);

    if (factor > max_inc)
      factor = max_inc;
    if (factor < max_dec)
      factor = max_dec;

    delta_t = factor * delta_t;


    if (delta_t < this->_min_time_step)
      throw SolveFailedException("Transient solver reached minimum time step");

    // Store the number of linear iterations required to
    // solve and the final residual.
    this->_n_linear_iterations = rval.first;
    this->_final_linear_residual = rval.second;

    // Update the system after the solve
    this->update();

    // this->time = _target_time;

    lin_solver->clear();
  }

  std::ostringstream os;
  os << "adaptive time stepping statistics:\n";
  os << "  current time step        : " << delta_t << "\n";
  os << "  relative time step error : " << rel_err << "\n";
  Messages::info(os.str());

  _time_step = delta_t;
}



template <class Base>
void
TiberTransientSystem<Base>::trapezoidal(void)
{

  // the full time step to achieve
  double delta_t = _target_time - this->time;

  // set the systems time to the target, for assembly
  //double old_time = this->time;
  this->time = _target_time;

  NumericVector<libMesh::Number>& A = this->get_vector("t_weight");

  // I think this should be always true, here
  if (this->assemble_before_solve)
  {
    A.zero();
    this->assemble();
  }

  // We have now all matrices/vectors at the new time step.
  // We use this as a temporary simplification. In principle,
  // for time-dependent parameters, or in nonlinear case,
  // they should refer to the two different times.

  // Get a reference to the EquationSystems
  const libMesh::EquationSystems& es =
    this->get_equation_systems();

  TiberLinearSolver* lin_solver =
    static_cast<TiberLinearSolver*>(this->linear_solver.get());

  if (get_options().has_submodel("linear_solver"))
  {
    auto lin_opts = get_options().submodels_begin("linear_solver");
    lin_solver->set_options(lin_opts->second);
  }

  double lin_rel_tol = lin_solver->get_linear_rtol();
  int lin_max_it = lin_solver->get_linear_max_it();

  //set the boundary values
  for (auto i = A.first_local_index(); i < A.last_local_index(); ++i)
  {
    if (A.el(i) == 0.0)
    {
      double val = this->rhs->el(i) / (*this->matrix)(i,i);
      this->solution->set(i, val);
    }
  }
  this->solution->close();

  // setup the final rhs
  // (2/delta_t*A - K)*c_old + 2*f
  //
  this->solution->scale(-1.0);
  this->rhs->scale(2.0);
  this->rhs->add_vector(*this->solution, *this->matrix);
  (*this->solution) *= A;
  this->rhs->add(-2.0/delta_t, *this->solution);

  // setup the final system matrix
  // 2/delta_t*A + K
  // A is the diagonal, so we have to add element wise
  //
  for (auto i = this->matrix->row_start(); i < this->matrix->row_stop(); ++i)
  {
    double val = 2.0/delta_t * A.el(i);
    this->matrix->add(i, i, val);
  }

  this->matrix->close();
  this->rhs->close();
  //this->matrix->print_matlab("K_solver.m");
  //this->rhs->print_matlab("F_solver.m");



  // Solve the linear system
  const std::pair<unsigned int, Real> rval;
  if (this->have_matrix("Preconditioner"))
    lin_solver->solve(*this->matrix, this->get_matrix("Preconditioner"),
        *this->solution, *this->rhs, lin_rel_tol, lin_max_it);
  else
    lin_solver->solve(*this->matrix, *this->solution, *this->rhs, lin_rel_tol, lin_max_it);


  // Store the number of linear iterations required to
  // solve and the final residual.
  this->_n_linear_iterations   = rval.first;
  this->_final_linear_residual = rval.second;

  // Update the system after the solve
  this->update();

  //this->time = _target_time;

  lin_solver->clear();
}


template <class Base>
void
TiberTransientSystem<Base>::user_initialization(void)
{
  ModelOptions solver_opts;
  auto lin_opts = get_options().submodels_begin("linear_solver");
  if (lin_opts != get_options().submodels_end("linear_solver"))
    solver_opts = lin_opts->second;
  
  TiberLinearSolver* lin_solver = TiberLinearSolver::create(this->comm(), solver_opts);
  libMesh::TransientSystem<Base>::linear_solver.reset();
  libMesh::TransientSystem<Base>::linear_solver = std::unique_ptr<libMesh::LinearSolver<Real> >(lin_solver);

  // we add a vector for the weight function of the time derivative
  this->add_vector("t_weight");
}




template <class Base>
std::string
TiberTransientSystem<Base>::system_type(void) const
{
  return "TiberTransient";
}


template <class Base>
libMesh::NumericVector<double>&
TiberTransientSystem<Base>::get_solution_vector(void)
{
  return *Base::current_local_solution;
}


template <class Base>
libMesh::NumericVector<double>&
TiberTransientSystem<Base>::get_local_solution_vector(void)
{
  return *Base::solution;
}





template class TiberTransientSystem<libMesh::LinearImplicitSystem>;
