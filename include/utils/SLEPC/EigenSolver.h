#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <vector>
#include <string>


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

  //!diagonalizes matrix
  static int eig_value_problem_general(const SLEPCoptions& opt) ;

  //!has to be called at the beginning of tibecad
  static void slepc_init(void);
  
  //!has to be called at the end of tibecad
  static void slepc_done(void);

  static int number_of_converged_eigenvalues();

  static double get_eigenvalue( int i);

  static void get_eigen_vector( int i, std::vector<Complex>& eigen_vector);

  static int prepare_slepc(void);

  static int clear_slepc(void);


 private:



  



};

#endif
