#include "TiberTransientSystem.h"
#include "TiberLinearSolver.h"

#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"


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
TiberTransientSystem<Base>::solve(void)
{

  // if the time has not been updated, we can return back
  if (_target_time <= this->time)
    return;

  if (this->assemble_before_solve)
    this->assemble();

  // Get a reference to the EquationSystems
  const libMesh::EquationSystems& es =
    this->get_equation_systems();

  TiberLinearSolver* lin_solver =
    static_cast<TiberLinearSolver*>(this->linear_solver.get());
  lin_solver->set_options(get_options());

  double lin_rel_tol = lin_solver->get_linear_rtol();
  int lin_max_it = lin_solver->get_linear_max_it();

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

  this->time = _target_time;

  lin_solver->clear();
}


template <class Base>
void
TiberTransientSystem<Base>::user_initialization(void)
{
  TiberLinearSolver* lin_solver = TiberLinearSolver::create(this->comm(), get_options());
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
