// $Id$

#ifndef _TIBERNONLINEARSYSTEM_H_
#define _TIBERNONLINEARSYSTEM_H_

#include "tibercad/solver/TiberEqSystem.h"
#include "tibercad/base/libMeshDefs.h"
#include "tibercad/base/tiber_dll.h"

#include "nonlinear_implicit_system.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"

//using namespace libMesh;


class XMonitor;


//! Base class for TiberCAD nonlinear systems
class TiberNonlinearSystem : public TiberEqSystem, public libMesh::NonlinearImplicitSystem
{

  public:

    typedef libMesh::NonlinearImplicitSystem Parent;

    //! The nonlinear solver implementations
    enum NonlinearSystemType
    {
      /*!
       * The TiberCAD internal Newton with line search
       */
      TIBER,

      /*! The PETSc nonlinear solver */
      PETSC,

      /*! The Bank and Rose nonlinear method (Bank, R. E. and Rose, D. J.:
       * "Global Approximate Newton Methods", Numerische Mathematik 37(2), 1981.
       */
      BANKROSE
    };

    
    
    //! The type of the assembly routine
    typedef void (*AssemblyRoutine)(const libMesh::NumericVector<Number> &X,
                                    libMesh::NumericVector<Number> *R,
                                    libMesh::SparseMatrix<Number> *J,
                                    libMesh::NonlinearImplicitSystem& system);
    
 

    //! Destructor
    virtual ~TiberNonlinearSystem(void) {};

   

    //! Create a nonlinear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param options the options for this system
     * \return a pointer to the newly created system
     */
    static TiberNonlinearSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, const ModelOptions& options);



    /*! \copydoc ImplicitSystem::clear() */
    virtual void clear(void);


    /*! \copydoc ImplicitSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void);


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const = 0;


    //! Get the solution vector including ghost values
    virtual libMesh::NumericVector<double>& get_solution_vector(void);

    //! Get the local solution vector (without ghost values)
    virtual libMesh::NumericVector<double>& get_local_solution_vector(void);


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
    TiberNonlinearSystem(libMesh::EquationSystems& es,
                         const std::string& name,
                         const unsigned int number);
   

    //! The assembly routine
    AssemblyRoutine _assemble;

    //! The real solve method, to be reimplemented in derived classes
    virtual void do_solve(void) = 0;

 
    //! The nonlinear iterations
    unsigned int _n_nonlin_iterations;


    //! The final residual norm
    double _final_residual_norm;


    //! The last Newton step size
    double _last_step_size;


    //! Get the pointer to the X monitor
    /*!
     * Use this with caution, especially don't forget to check for NULL pointer
     */
    XMonitor* get_xmonitor(void);


    //! Add a point to the X monitor
    /*!
     * \param iteration the iteration number
     * \param err the error
     * \param logarithm if \c true, plot \c log10(error)
     */
    void draw_point(double iteration, double error, bool logarithm = true);


  private:

    //! The X monitor
    XMonitor* _xmonitor;

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
  ImplicitSystem::clear();
}



inline
void
TiberNonlinearSystem::reinit(void)
{
  ImplicitSystem::reinit();
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
XMonitor*
TiberNonlinearSystem::get_xmonitor(void)
{
  return _xmonitor;
}


#endif // _TIBERNONLINEARSYSTEM_H_
