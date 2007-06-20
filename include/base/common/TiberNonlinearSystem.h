// $Id$

#ifndef _TIBERNONLINEARSYSTEM_H_
#define _TIBERNONLINEARSYSTEM_H_

//#include "ModelOptions.h"

#include "implicit_system.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"


template<typename> class LinearSolver;
class EquationSystems;


//! A generic class to solve nonlinear systems
/*!
 * This implementation uses an approximate Newton scheme to solve a
 * nonlinear system
 * \f[g'(x^k)\delta x^k = -g(x^k)\f]
 * with update strategy
 * \f[x^{k+1} = x^k + t_k\delta x^k
 * 
 * The relaxation factors \f$t_k\f$ are currently calculated using a
 * standard line search algorithm.
 */
class TiberNonlinearSystem : public ImplicitSystem
{

  public:

    //! The nonlinear solver implementations
    enum NonlinearSystemType
    {
      /*!
       * The TiberCAD internal Newton with line search
       */
      TIBER,

      /*! The PETSc nonlinear solver */
      PETSC
    };

    
    //! The type of norms
    enum NormType
    {
      MAX_NORM,  //< the maximum norm
      l2_NORM,   //< the l2 norm
    };

    
    //! The type of the assembly routine
    typedef void (*AssemblyRoutine)(const NumericVector<Number> &X,
                                    NumericVector<Number> *R,
                                    SparseMatrix<Number> *J);
    
 

    //! Destructor
    virtual ~TiberNonlinearSystem(void) {};


    //! To create a nonlinear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param type the type of system to create
     * \return a reference to the newly created system
     */
    static TiberNonlinearSystem& create_nonlinear_system(EquationSystems& es,
        const std::string& sysname, NonlinearSystemType type,
        const std::string& linear_solver = "petsc");

    
    //! To create a nonlinear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param type the type of system to create as string
     * \return a reference to the newly created system
     */
    static TiberNonlinearSystem& create_nonlinear_system(EquationSystems& es,
        const std::string& sysname, const std::string& type,
        const std::string& linear_solver = "petsc");


    /*! \copydoc ImplicitSystem::clear() */
    virtual void clear(void);


    /*! \copydoc ImplicitSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void) = 0;


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const;


    //! Get the solution vector
    virtual NumericVector<double>& get_solution_vector(void) = 0;


    //! Set the linear solver parameters
    void set_linear_solver_params(double relative_tolerance,
                                  double absolute_tolerance,
                                  unsigned int maximum_iterations);


    //! Set the linear solver type
    void set_linear_solver_type(SolverType solver_type, PreconditionerType pc);


    //! Set the nonlinear solver parameters
    void set_nonlinear_solver_params(double relative_tolerance,
                                     double absolute_tolerance,
                                     double step_tolerance,
                                     double max_step,
                                     unsigned int maximum_iterations);


    //! Attach the assembly routine
    void attach_assembly_routine(AssemblyRoutine assembly);


    //! Get the number of nonlinear iterations
    unsigned int n_nonlinear_iterations(void) const;


    //! Get the norm of the final residual
    double final_residual_norm(void) const;


    //! Get the norm of the last approximate Newton step
    double last_step_size(void) const;



  protected:

    //! Constructor
    TiberNonlinearSystem(EquationSystems& es,
                         const std::string& name,
                         const unsigned int number);
   

    //! The assembly routine
    AssemblyRoutine _assemble;

 
    //! The nonlinear iterations
    unsigned int _n_nonlin_iterations;


    //! The final residual norm
    double _final_residual_norm;


    //! The last Newton step size
    double _last_step_size;


    //! The tolerance for the nonlinear step size
    double _nonlin_step_tol;


    //! The relative tolerance for the residual norm
    double _nonlin_rel_tol;


    //! The absolute tolerance for the residual norm
    double _nonlin_abs_tol;


    //! The nonlinear maximum number of iterations
    unsigned int _nonlin_max_it;


    //! The linear relative tolerance
    double _lin_tol;


    //! The linear absolute tolerance
    double _lin_abs_tol;


    //! The linear solver maximum number of iterations
    unsigned int _lin_max_it;


    //! The maximum step size
    double _max_step_size;


    //! The linear solver type
    SolverType _solver_type;


    //! The preconditioner type
    PreconditionerType _preconditioner_type;
   

    //! The linear solver
    std::string _linear_solver;


  private:

    //! The parent class type
    typedef ImplicitSystem Parent;

};



//
// inline methods
//


inline
void
TiberNonlinearSystem::attach_assembly_routine(AssemblyRoutine assembly)
{
  _assemble = assembly;
}


inline
void
TiberNonlinearSystem::clear(void)
{
  Parent::clear();
}



inline
void
TiberNonlinearSystem::reinit(void)
{
  Parent::reinit();
}



    
inline
unsigned int
TiberNonlinearSystem::n_nonlinear_iterations(void) const
{
  return _n_nonlin_iterations;
}


inline
double
TiberNonlinearSystem::final_residual_norm(void) const
{
  return _final_residual_norm;
}


inline
double
TiberNonlinearSystem::last_step_size(void) const
{
  return _last_step_size;
}


inline
std::string
TiberNonlinearSystem::system_type(void) const
{
  return "TiberNonlinear";
}


inline
void
TiberNonlinearSystem::set_linear_solver_params(double relative_tolerance,
    double absolute_tolerance, unsigned int maximum_iterations)
{
  _lin_tol = relative_tolerance;
  _lin_abs_tol = absolute_tolerance;
  _lin_max_it = maximum_iterations;
}


inline
void
TiberNonlinearSystem::set_nonlinear_solver_params(double relative_tolerance,
    double absolute_tolerance, double step_tolerance, double max_step,
    unsigned int maximum_iterations)
{
  _nonlin_rel_tol = relative_tolerance;
  _nonlin_abs_tol = absolute_tolerance;
  _nonlin_step_tol = step_tolerance;
  _max_step_size = max_step;
  _nonlin_max_it = maximum_iterations;
}


inline
void
TiberNonlinearSystem::set_linear_solver_type(SolverType solver_type,
    PreconditionerType pc)
{
  _solver_type = solver_type;
  _preconditioner_type = pc;
}


#endif // _TIBERNONLINEARSYSTEM_H_
