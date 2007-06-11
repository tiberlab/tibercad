// $Id$

#include "PardisoLinearSolver.h"


#include "petsc_vector.h"
#include "petsc_matrix.h"



#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petscksp.h>
}
#else
# include <petscversion.h>
# include <petscksp.h>
#endif



template <typename T>
void
PardisoLinearSolver<T>::clear(void)
{
}


template <typename T>
void
PardisoLinearSolver<T>::init(void)
{
}


template <typename T>
std::pair<unsigned int, Real> 
PardisoLinearSolver<T>::solve(SparseMatrix<T>&  matrix_in,
			     SparseMatrix<T>&  precond_in,
			     NumericVector<T>& solution_in,
			     NumericVector<T>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{
}


// Explicit instantiation
template class PardisoLinearSolver<Number>;
 

