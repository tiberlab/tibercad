/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file PardisoLinearSolver.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef TC_PARDISOLINEARSOLVER_H
#define TC_PARDISOLINEARSOLVER_H


#include "TiberLinearSolver.h"
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
    PardisoLinearSolver(const ModelOptions& options);


    //! Destructor.
    virtual ~PardisoLinearSolver(void);


    //! Release all memory and clear data structures.
    virtual void clear(void);

    //! Initialize data structures.
    virtual void init(void);



  protected:

    //! Call the linear solver specifying explicitly the preconditioner matrix
    virtual std::pair<unsigned int, Real> 
      do_solve (SparseMatrix<Number>  &matrix,
          SparseMatrix<Number>  &preconditioner,
          NumericVector<Number> &solution,
          NumericVector<Number> &rhs);



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
    static void _checkerr(int errorcode);

    //! Pardiso solver interface
    void solve_pardiso(double *mat, int *ia, int *ja, double *b, double *x, int n);

   
};



//
// inline members
//


inline
void
PardisoLinearSolver::_checkerr(int errorcode)
{
  if (errorcode != 0)
    throw PardisoSolverException(errorcode);
}





inline
PardisoLinearSolver::PardisoLinearSolver(const ModelOptions& options) :
  TiberLinearSolver(options)
{
}



inline
PardisoLinearSolver::~PardisoLinearSolver(void)
{
  this->clear();
}





#endif // TC_PARDISOLINEARSOLVER_H
