// $Id$


#ifndef _TIBERLINEARSOLVER_H_
#define _TIBERLINEARSOLVER_H_

#include "TypeDefs.h"

// Libmesh includes
#include "linear_solver.h"

class ModelOptions;


//! The TiberCAD linear solver interface
/*!
 * This class provides the TiberCAD interface to
 * linear solvers.
 * It is derived from the libmesh LinearSolver class
 *
 */
class TiberLinearSolver : public LinearSolver<Number>
{

  public:

    //!  Constructor. Initializes  data structures
    TiberLinearSolver(void);


    //! Destructor.
    virtual ~TiberLinearSolver(void);

    
    //! Create a linear solver
    static TiberLinearSolver* create(const std::string& type);

    
    //! Create a linear solver
    static TiberLinearSolver* create(const ModelOptions& options);
    

    //! Release all memory and clear data structures.
    virtual void clear(void) { };


    //! Initialize data structures if not done so already.
    virtual void init(void) = 0;


    //! Set the options for the linear solver
    /*!
     * \deprecated This method will disappear
     */
    void set_ksp_options(double rtol = 1e-6, unsigned int max_it = 1000);


    //! Set the options for the linear solver
    /*!
     * \deprecated This method will disappear
     */
    void set_ksp_options(double rtol, double atol, unsigned int max_it = 1000);


    //! Set options
    /*!
     * Call this method before calling solve().
     *
     * Unspecified options are set to there default values,
     * \em not to their current values!
     */
    void set_options(const ModelOptions& options);


    //! Call the  solver.
    /*!
     * It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     *
     * This method is used only for compatibility with libMesh.
     */    
    virtual std::pair<unsigned int, Real> 
      solve(SparseMatrix<Number> &matrix,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs,
          const double tol,
          const unsigned int m_its)
      {
        ignore_unused_variable(tol);
        ignore_unused_variable(m_its);
        return this->solve(matrix, matrix, solution, rhs);
      }


    //! Call the linear solver specifying explicitly the preconditioner matrix
    /*!
     * This method is used only for compatibility with libMesh.
     */
    virtual std::pair<unsigned int, Real> 
      solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs,
          const double tol,
          const unsigned int m_its)
      {
        ignore_unused_variable(tol);
        ignore_unused_variable(m_its);
        return this->solve(matrix, preconditioner, solution, rhs);
      }


    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    virtual std::pair<unsigned int, Real>
      solve(SparseMatrix<Number> &matrix,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs)
      {
        return this->solve(matrix, matrix, solution, rhs);
      }



    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    virtual std::pair<unsigned int, Real>
      solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs) = 0;


    //! Get the relative linear tolerance
    double get_linear_rtol(void) const;

    //! Get the absolute linear tolerance
    double get_linear_atol(void) const;

    //! Get the maximum number of iterations
    int get_linear_max_it(void) const;

    //! Set the relative linear tolerance
    void set_linear_rtol(double rtol);

    //! Set the absolute linear tolerance
    void set_linear_atol(double atol);

    //! Set the maximum number of iterations
    void set_linear_max_it(int max_it);


    //! Get simulation name
    const std::string& get_simulation_name(void) const;



  protected:

    //! Parse the options for solver specific stuff
    virtual void parse_options(const ModelOptions& options);



  private:

    //! The relative linear tolerance
    double _linear_rtol;

    //! The absolute linear tolerance
    double _linear_atol;
      
    //! The maximum number of iterations
    int _linear_max_it;

    //! The name of the associated simulation
    std::string _sim_name;

};



//
// inline members
//



inline
TiberLinearSolver::~TiberLinearSolver(void)
{
}


///*
inline
void
TiberLinearSolver::set_ksp_options(double rtol,
    unsigned int max_it)
{
  set_ksp_options(rtol, 1e-50, max_it);
}



inline
void
TiberLinearSolver::set_ksp_options(double rtol, double atol,
    unsigned int max_it)
{
  _linear_rtol = rtol;
  _linear_atol = atol;
  _linear_max_it = max_it;
}
//*/

inline
double
TiberLinearSolver::get_linear_rtol(void) const
{
  return _linear_rtol;
}


inline
double
TiberLinearSolver::get_linear_atol(void) const
{
  return _linear_atol;
}


inline
int
TiberLinearSolver::get_linear_max_it(void) const
{
  return _linear_max_it;
}


inline
void
TiberLinearSolver::set_linear_rtol(double rtol)
{
  _linear_rtol = rtol;
}


inline
void
TiberLinearSolver::set_linear_atol(double atol)
{
  _linear_atol = atol;
}


inline
void
TiberLinearSolver::set_linear_max_it(int max_it)
{
  _linear_max_it = max_it;
}


inline
void
TiberLinearSolver::parse_options(const ModelOptions& options)
{
  static_cast<const void*>(&options);
}

inline
const std::string&
TiberLinearSolver::get_simulation_name(void) const
{
  return _sim_name;
}



#endif // _TIBERLINEARSOLVER_H_
