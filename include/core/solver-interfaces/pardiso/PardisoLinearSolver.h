// $Id$


#ifndef _PARDISOLINEARSOLVER_H_
#define _PARDISOLINEARSOLVER_H_

#include "LinearSolverException.h"


// Libmesh includes
#include "linear_solver.h"



//! The TiberCAD linear solver interface to PARDISO
/*!
 * This class provides the TiberCAD interface to the PARDISO
 * linear solver.
 * It is derived from the libmesh LinearSolver class
 *
 */
template <typename T>
class PardisoLinearSolver : public LinearSolver<T>
{

  public:

    //!  Constructor. Initializes Pardiso data structures
    PardisoLinearSolver(void);


    //! Destructor.
    virtual ~PardisoLinearSolver(void);


    //! Release all memory and clear data structures.
    virtual void clear(void);


    //! Initialize data structures if not done so already.
    virtual void init(void);


    //! Call the Pardiso solver.
    /*!
     * It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */    
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<T>  &matrix_in,
          NumericVector<T> &solution_in,
          NumericVector<T> &rhs_in,
          const double tol,
          const unsigned int m_its)
      {
        return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
      }


    //! Call the linear solver specifying explicitly the preconditioner matrix
    virtual std::pair<unsigned int, Real> 
      solve (SparseMatrix<T>  &matrix,
          SparseMatrix<T>  &preconditioner,
          NumericVector<T> &solution,
          NumericVector<T> &rhs,
          const double tol,
          const unsigned int m_its);



  private:

        //! Maximum number of numerical factorizations. 
	int maxfct; 

	//! Which factorization to use. 
        int  mnum; 

	//! Initialize error flag 
        int  error;  

	//! Print statistical information in file */
        int  msglvl; 

        //! Real unsymmetric matrix 
	int mtype; 

        //! Number of right hand sides.
  	int nrhs; 

        //! Parameters of Pardiso solver
	int iparm[64]; 





    //! Check PETSc error code
    static void _checkerr(int errorcode) throw (LinearSolverException);

};



//
// inline members
//


template <typename T>
inline
void
PardisoLinearSolver<T>::_checkerr(int errorcode) throw (LinearSolverException)
{
  if (errorcode != 0)
    throw(LinearSolverException("Petsc error in PardisoLinearSolver"));
}





template <typename T>
inline
PardisoLinearSolver<T>::PardisoLinearSolver(void)
{
}



template <typename T>
inline
PardisoLinearSolver<T>::~PardisoLinearSolver(void)
{
  this->clear();
}





#endif // _PARDISOLINEARSOLVER_H_
