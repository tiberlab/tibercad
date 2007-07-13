// $Id$


#ifndef _PARDISOLINEARSOLVER_H_
#define _PARDISOLINEARSOLVER_H_


#include "TiberLinearSolver.h"
#include "LinearSolverException.h"
#include "PardisoSolverException.h"



//! The TiberCAD linear solver interface to PARDISO
/*!
 * This class provides the TiberCAD interface to the PARDISO
 * linear solver.
 * It is derived from the libmesh LinearSolver class
 *
 */
class PardisoLinearSolver : public TiberLinearSolver
{

  public:

    //!  Constructor. Initializes Pardiso data structures
    PardisoLinearSolver(void);


    //! Destructor.
    virtual ~PardisoLinearSolver(void);


    //! Release all memory and clear data structures.
    virtual void clear(void);

    //! Initialize data structures.
    virtual void init(void);


    //! Call the linear solver specifying explicitly the preconditioner matrix
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<Number>  &matrix,
          SparseMatrix<Number>  &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs,
          const double tol,
          const unsigned int m_its);



  private:

    //! Maximum number of numerical factorizations. 
    int maxfct; 

    //! Which factorization to use. 
    int  mnum; 

    //! Initialize error flag 
    int  error;  

    //! Print statistical information in file
    int  msglvl; 

    //! Real unsymmetric matrix 
    int mtype; 

    //! Number of right hand sides.
    int nrhs; 

    //! Parameters of Pardiso solver
    int iparm[64]; 

    //! Check PETSc error code
    static void _checkerr(int errorcode) throw (LinearSolverException);

    //! Pardiso solver interface
    void solve_pardiso(double *mat, int *ia, int *ja, double *b, double *x, int n);

   
};



//
// inline members
//


inline
void
PardisoLinearSolver::_checkerr(int errorcode) throw (LinearSolverException)
{
  if (errorcode != 0)
  {
    std::ostringstream s;
    s << "Petsc error in PardisoLinearSolver: " << errorcode;
    throw(LinearSolverException(s.str()));
  }
}





inline
PardisoLinearSolver::PardisoLinearSolver(void)
{
}



inline
PardisoLinearSolver::~PardisoLinearSolver(void)
{
  this->clear();
}





#endif // _PARDISOLINEARSOLVER_H_
