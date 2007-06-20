// $Id$


#ifndef _TIBERLINEARSOLVER_H_
#define _TIBERLINEARSOLVER_H_


// Libmesh includes
#include "linear_solver.h"



//! The TiberCAD linear solver interface to PETSc
/*!
 * This class provides the TiberCAD interface to PETSc
 * iterative solvers.
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
    

    //! Release all memory and clear data structures.
    virtual void clear(void) { };


    //! Initialize data structures if not done so already.
    virtual void init(void) = 0;


    //! Set the options for the linear solver
    void set_ksp_options(double rtol = 1e-6, unsigned int max_it = 1000);


    //! Set the options for the linear solver
    void set_ksp_options(double rtol, double atol, unsigned int max_it = 1000);


    //! Call the  solver.
    /*!
     * It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */    
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<Number>  &matrix_in,
          NumericVector<Number> &solution_in,
          NumericVector<Number> &rhs_in,
          const double tol,
          const unsigned int m_its)
      {
        return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
      }


    //! Call the linear solver specifying explicitly the preconditioner matrix
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<Number>  &matrix,
          SparseMatrix<Number>  &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs,
          const double tol,
          const unsigned int m_its) = 0;


  protected:

    double get_linear_rtol();
    double get_linear_atol();
    int get_linear_max_it();

  private:

    double _linear_rtol;
    double _linear_atol;
    int _linear_max_it;

};



//
// inline members
//


inline
TiberLinearSolver::TiberLinearSolver(void)
  : _linear_rtol(1e-6),
    _linear_atol(1e-50),
    _linear_max_it(500)
{
}



inline
TiberLinearSolver::~TiberLinearSolver(void)
{
}



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


inline
double
TiberLinearSolver::get_linear_rtol()
{
  return _linear_rtol;
}


inline
double
TiberLinearSolver::get_linear_atol()
{
  return _linear_atol;
}


inline
int
TiberLinearSolver::get_linear_max_it()
{
  return _linear_max_it;
}



#endif // _TIBERLINEARSOLVER_H_
