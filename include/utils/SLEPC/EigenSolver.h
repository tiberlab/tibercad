#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <vector>
#include <string>
#include <complex>

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
    
    std::string st_ksp_type; //!< Liner system solution method 
    
    double spectrum_shift;

    bool read_matrix_from_file;

    bool matrix_output;

    std::string pc_type; //!<preconditioner name


    bool monitor;   //<! activates monitor if true
    
  };

  //!solves general eigenvalue problem (Hx = gSx) matrix
  static int eig_value_problem_general(const SLEPCoptions& opt) ;


  //!solves  eigenvalue problem (Hx = gx)  matrix
  static int eig_value_problem(const SLEPCoptions& opt) ;

  //!has to be called at the beginning of tibecad
  static void slepc_init(int argc1, char** argv1);
  
  //!has to be called at the end of tibecad
  static void slepc_done(void);

  //!returns number of converged eigenvalues
  static int number_of_converged_eigenvalues();

  //!returns  eigenvalue number i (starting from 0) 
  static double get_eigenvalue( int i);

  //!returns  eigenvector number i (starting from 0) 
  static void get_eigen_vector( int i, std::vector<Complex>& eigen_vector);

  //!has to be called before solving the eigenvalue problem
  static int prepare_slepc(void);

  //!has to be called after solving the eigenvalue problem
  static int clear_slepc(void);

  //!init H_matrix
  /*!
    \param size  matrix size
   */
  static int init_H_matrix(unsigned int size);

  //!init S_matrix
  /*!
    \param size  matrix size
   */
  static int init_S_matrix(unsigned int size);


  //!allocates memory for H matrix (non-parallel version!!!) 
  /*!
    \param matrix_size size of the matrix 
    \param non_zeros numbers of non-zero columns in each raw
  */
  static int preallocate_H_matrix(unsigned int matrix_size,  int*  non_zeros);

  //!allocates memory for S matrix (non-parallel version!!!) 
  /*!
    \param matrix_size size of the matrix 
    \param non_zeros numbers of non-zero columns in each raw
  */
  static int preallocate_S_matrix(unsigned int matrix_size,  int*  non_zeros);
  

  //closes H matrix
  static int finalize_H_assembly(void);

  //close S matrix
  static int finalize_S_assembly(void);

  //!pass a row into H matrix
  /*!
    \param row row number
    \param columns numbers of columns
    \param value values of matrix elements
  */
  static int insert_H_row( int row, const std::vector<unsigned int>& columns, const std::vector<Complex>& value);

  //!pass a row into S matrix
  /*!
    \param row row number
    \param columns numbers of columns
    \param value values of matrix elements
  */
  static int insert_S_row( int row, const std::vector<unsigned int>& columns, const std::vector<Complex>& value);



  //!set initial vector
  /*!
    \param  initial_vector initial vector 
  */
  static void set_initial_vector( const std::vector<Complex>& initial_vector);


 private:

  static int do_solve(const SLEPCoptions& opt);

  
  static  int _size_of_matrix;


};

#endif
