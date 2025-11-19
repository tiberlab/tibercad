#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include "tibercad/base/tiber_dll.h"

#include <vector>
#include <string>
#include <complex>
#include <mpi.h>

namespace libMesh {
  typedef std::complex<double> Complex;
}
using libMesh::Complex;


//! SLEPc interface class 
class TBDLEXPORT EigenSolver
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
    
    Complex spectrum_shift;

    bool read_matrix_from_file;

    bool matrix_output;

    bool use_deflation_space;

    std::string pc_type; //<! preconditioner name

    std::string solver_package{"mumps"}; //<! solver package, e.g. mumps, petsc

    std::string spectral_trans;

    bool monitor;   //<! activates monitor if true

    double spectrum_inversion_tolerance; //<! toterance for spectrum inversion

  };

  //! Eigenvalue problem types
  enum EVPType
  {
    NORMAL = 0, //! normal EVP Ax = kx
    GENERALIZED = 1 //! generalized EVP Ax = kBx
  };



  //! solves general eigenvalue problem (Hx = gSx) matrix
  static int eig_value_problem_general(const SLEPCoptions& opt);

  //! solve a normal or general eigenvalue problem
  /*!
   * normal eigenvalue problem is Hx = gx
   * general eigenvalue problem is Hx = gSx
   *
   * \param opt the solver options
   * \param evp_type eigenproblem type, \c NORMAL or \c GENERALIZED
   */
  static int eig_value_problem(const SLEPCoptions& opt,
                               EVPType evp_type = NORMAL);

  //!has to be called at the beginning of tibecad
  static void slepc_init(int argc1, char** argv1, MPI_Comm comm);

  static void set_deflation_space(
      const std::vector<const std::vector<Complex>*>& solutions);

  //!has to be called at the end of tibecad
  static void slepc_done(void);

  //!returns number of converged eigenvalues
  static int number_of_converged_eigenvalues();

  //!returns  eigenvalue number i (starting from 0)
  static double get_eigenvalue( int i);

  //!returns  eigenvector number i (starting from 0)
  static void get_eigen_vector( int i, std::vector<Complex>& eigen_vector);

  //!has to be called before solving the eigenvalue problem
  static int prepare_slepc(MPI_Comm comm);

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


  //! Allocates memory for H or S matrix
  /*!
    \param matrix the matrix type, 'H' or 'S'
    \param N the global size
    \param n the local size
    \param d_nnz the nonzeros on each local row in the diagonal block 
    \param o_nnz the nonzeros on each local row in the off-diagonal block 
  */
  static int preallocate_matrix(const char matrix, unsigned int N, unsigned int n,
      std::vector<int>& d_nnz, std::vector<int>& o_nnz);


  //! Finalize assmebly of H or S matrix
  static void finalize_matrix_assembly(const char matrix);

  //! pass a row into H or S matrix
  /*!
    \param matrix matrix type. 'H' or 'S'
    \param row row number
    \param columns numbers of columns
    \param value values of matrix elements
    \param indexing_base the base for indexing in \c columns
  */
  static void insert_matrix_row(const char matrix, int row,
      const std::vector<unsigned int>& columns,
      const std::vector<Complex>& value, int indexing_base = 0);

  static void print_options(const SLEPCoptions& opt);

  //!set initial vector
  /*!
    \param  initial_vector initial vector
  */
  static void set_initial_vector( const std::vector<Complex>& initial_vector);

  //! get actual shift used in calculation
  static double get_shift(void);

  //! Check whether matrices are Hermitian
  static bool check_matrices(double tol, bool verbose);
  static bool check_matrices(void);


 private:

  static int do_solve(const SLEPCoptions& opt);


  static  int _size_of_matrix;


};



#endif

