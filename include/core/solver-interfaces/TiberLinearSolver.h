// $Id$


#ifndef _TIBERLINEARSOLVER_H_
#define _TIBERLINEARSOLVER_H_

//#include "TypeDefs.h"
#include "TiberModelObject.h"

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
class TiberLinearSolver : public TiberModelObject, public LinearSolver<Number>
{

  public:

    //!  Constructor. Initializes  data structures
    TiberLinearSolver(const ModelOptions& options);


    //! Destructor.
    virtual ~TiberLinearSolver(void);

    
    //! Create a linear solver
    static TiberLinearSolver* create(const ModelOptions& options);
    

    //! Release all memory and clear data structures.
    virtual void clear(void) { };


    //! Initialize data structures if not done so already.
    virtual void init(void) = 0;



    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    std::pair<unsigned int, Real>
      solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs)
      {
        parse_options();
        return this->do_solve(matrix, matrix, solution, rhs);
      }


    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    std::pair<unsigned int, Real>
      solve(SparseMatrix<Number> &matrix,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs)
      {
        parse_options();
        return this->do_solve(matrix, matrix, solution, rhs);
      }


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
          const double,
          const unsigned int)
      {
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
          const double,
          const unsigned int)
      {
        return this->solve(matrix, preconditioner, solution, rhs);
      }



    //! This is not supported in TiberCAD
    virtual std::pair<unsigned int, Real>
      solve(const ShellMatrix<Number>&,
        NumericVector<Number>&,
        NumericVector<Number>&,
          const double,
          const unsigned int);
  


    //! This is not supported in TiberCAD
    virtual std::pair<unsigned int, Real>
      solve(const ShellMatrix<Number>&,
          const SparseMatrix<Number>&,
          NumericVector<Number>&,
          NumericVector<Number>&,
          const double,
          const unsigned int);
  



    //! Print a message on convergence
    virtual void print_converged_reason() {};


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
    std::string get_simulation_name(void) const;



  protected:

    //! Parse the options for solver specific stuff
    /*!
     * Calls do_parse_options()
     */
    void parse_options(void);


    //! Parse the options for solver specific stuff
    virtual void do_parse_options(void);


    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    virtual std::pair<unsigned int, Real>
      do_solve(SparseMatrix<Number> &matrix,
          SparseMatrix<Number> &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs) = 0;


    //! Solve the linear system
    /*!
     * In TiberCAD, this method should be used instead of the
     * ones defined in libMesh
     */
    virtual std::pair<unsigned int, Real>
      do_solve(SparseMatrix<Number> &matrix,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs)
      {
        return this->do_solve(matrix, matrix, solution, rhs);
      }


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
  get_options().set_option("relative_tolerance", rtol);
  _linear_rtol = rtol;
}


inline
void
TiberLinearSolver::set_linear_atol(double atol)
{
  get_options().set_option("absolute_tolerance", atol);
  _linear_atol = atol;
}


inline
void
TiberLinearSolver::set_linear_max_it(int max_it)
{
  get_options().set_option("max_iterations", max_it);
  _linear_max_it = max_it;
}



inline
void
TiberLinearSolver::do_parse_options(void)
{
}


#endif // _TIBERLINEARSOLVER_H_
