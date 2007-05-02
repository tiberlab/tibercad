#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <vector>
#include <string>
#include "slepceps.h"

typedef std::complex<double> Complex;

//! SLEPc interface class 
class EigenSolver
{
 public:
  //!parameters for solver
  struct SLEPCoptions
  {
    unsigned int ev_number;

    std::string H_file_name;
    
    std::string S_file_name;
    
    std::string solver_type;
    
    double eps_tolerance;
    
    unsigned int eps_max_it;
    
    double st_ksp_rtol;
    
    std::string st_ksp_type;
    
    double spectrum_shift;
    
  };


  static int eig_value_problem_general(const SLEPCoptions& opt) ;


  static void slepc_init(void);

  static void slepc_done(void);

  static inline int number_of_converged_eigenvalues();

  static inline double get_eigenvalue( int i);

  static inline void get_eigen_vector( int i, std::vector<Complex>& eigen_vector);


 private:
  //!Matrix to be diagonalized
  static Mat         A;

  //!S matrix
  static Mat         B;  

  //! eigenproblem solver context 
  static EPS         eps;             
  

  



};
//--------------------------------------------------------------//
inline  int EigenSolver::number_of_converged_eigenvalues()
{
  int ierr, nconv; 
  ierr =  EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  return(nconv);
}
//--------------------------------------------------------------//
inline double EigenSolver::get_eigenvalue( int i)
{
  int ierr;
  PetscScalar ev, ev_i;
/*  
  ierr = EPSGetValue(eps, i, &ev,  &ev_i);

  double eigen_value = PetscRealPart(ev);

  return(eigen_value);
 */
}
//----------------------------------------------------------------//


inline void EigenSolver::get_eigen_vector( int i, std::vector<Complex>& eigen_vector_out)
{
  int ierr, vec_size;
  PetscScalar kr, ki;
  Vec eigen_vector;

  EPSGetEigenpair(eps,i,&kr,&ki,eigen_vector,PETSC_NULL);
 
  VecGetSize(eigen_vector, &vec_size);
  
  eigen_vector_out.resize(vec_size);
  for (int j= 0; j < vec_size; j++)
  {
    //  eigen_vector_out[j] = Complex(  PetscRealPart(eigen_vector[j]), PetscImaginaryPart(eigen_vector[j]) );  
  }

}
#endif
