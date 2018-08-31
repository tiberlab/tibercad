// $Id$

#include "TiberMath.h"
#include "Messages.h"

#include <cmath>


using namespace std;

pair<double, double>
Distributions::fermi_dirac(double E, double kT)
{
  double f = 0, deriv = 0;
  double g = 1;
  double arg = -E / kT;
  if (arg > 50)
  {
    f = exp(-arg) / g;
    deriv = -f;
  }
  else if (arg < -50)
  {
    deriv = -g * exp(arg);
    f = 1 + deriv;
  }
  else
  {
    double expfac = g * exp(arg);
    double denom = 1.0 + expfac;
    f = 1.0 / denom;
    deriv = -expfac * f / denom;
  }

  return make_pair(f, -deriv / kT);
}

// added to calculate the second derivative
void
Distributions::fermi_dirac(std::vector<double>& result,
    double E, double kT)
{
  double f = 0, deriv = 0, deriv2 = 0;
  double g = 1;
  double arg = -E / kT;
  if (arg > 50)
  {

    //f = exp(-arg) / g;
    result[0] = exp(-arg) / g;
    if (result.size() > 1)
      result[1] = -result[0];
      //deriv = -f;
    if (result.size() > 2)
      result[2] = result[0];
      //deriv2 = f;
  }
  else if (arg < -50)
  {

    result[0] = 1 + -g * exp(arg);
    if (result.size() > 1)
      result[1] =-g * exp(arg);
    if (result.size() > 2)
      result[2] = - result[1];

  }
  else
  {

    double expfac = g * exp(arg);
    double denom = 1.0 + expfac;

    result[0] = 1.0 / denom;
    if (result.size() > 1)
      result[1] = -expfac * result[0] / denom;
    if (result.size() > 2)
      result[2] = result[1] + 2 * expfac * result[0] / denom / denom;


  }

  result[1] *= 1 / kT;
  result[2] *= 1 / kT / kT;

}


/*
#include "dense_matrix.h"
#include "dense_vector.h"

#include "petsc_macro.h"
EXTERN_C_FOR_PETSC_BEGIN
#include <petscblaslapack.h>
EXTERN_C_FOR_PETSC_END


void
TiberMath::svd(DenseMatrix<double>& matrix, DenseVector<double>& sigma)
{
  char JOBU = 'N';
  char JOBVT = 'N';
  std::vector<double> sigma_val;
  std::vector<double> U_val;
  std::vector<double> VT_val;

  //    M       (input) int*
  //            The number of rows of the matrix A.  M >= 0.
  // In C/C++, pass the number of *cols* of A
  int M = matrix.n();

  //    N       (input) int*
  //            The number of columns of the matrix A.  N >= 0.
  // In C/C++, pass the number of *rows* of A
  int N = matrix.m();

  int min_MN = (M < N) ? M : N;
  int max_MN = (M > N) ? M : N;

  //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //          On entry, the M-by-N matrix A.
  //          On exit,
  //          if JOBU = 'O',  A is overwritten with the first min(m,n)
  //                          columns of U (the left singular vectors,
  //                          stored columnwise);
  //          if JOBVT = 'O', A is overwritten with the first min(m,n)
  //                          rows of V**T (the right singular vectors,
  //                          stored rowwise);
  //          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
  //                          are destroyed.
  // Here, we pass &(_val[0]).

  //    LDA     (input) int*
  //            The leading dimension of the array A.  LDA >= max(1,M).
  int LDA = M;

  //  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //          The singular values of A, sorted so that S(i) >= S(i+1).
  sigma_val.resize( min_MN );

  //  LDU     (input) INTEGER
  //          The leading dimension of the array U.  LDU >= 1; if
  //          JOBU = 'S' or 'A', LDU >= M.
  int LDU = M;

  //  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
  //          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
  //          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
  //          if JOBU = 'S', U contains the first min(m,n) columns of U
  //          (the left singular vectors, stored columnwise);
  //          if JOBU = 'N' or 'O', U is not referenced.
  U_val.resize( LDU*M );

  //  LDVT    (input) INTEGER
  //          The leading dimension of the array VT.  LDVT >= 1; if
  //          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
  int LDVT = N;

  //  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
  //          V**T;
  //          if JOBVT = 'S', VT contains the first min(m,n) rows of
  //          V**T (the right singular vectors, stored rowwise);
  //          if JOBVT = 'N' or 'O', VT is not referenced.
  VT_val.resize( LDVT*N );

  //  LWORK   (input) INTEGER
  //          The dimension of the array WORK.
  //          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
  //          For good performance, LWORK should generally be larger.
  //
  //          If LWORK = -1, then a workspace query is assumed; the routine
  //          only calculates the optimal size of the WORK array, returns
  //          this value as the first entry of the WORK array, and no error
  //          message related to LWORK is issued by XERBLA.
  int larger = (3*min_MN+max_MN > 5*min_MN) ? 3*min_MN+max_MN : 5*min_MN;
  int LWORK  = (larger > 1) ? larger : 1;


  //  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
  //          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
  //          superdiagonal elements of an upper bidiagonal matrix B
  //          whose diagonal is in S (not necessarily sorted). B
  //          satisfies A = U * B * VT, so it has the same singular values
  //          as A, and singular vectors related by U and VT.
  std::vector<double> WORK( LWORK );

  //  INFO    (output) INTEGER
  //          = 0:  successful exit.
  //          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //          > 0:  if DBDSQR did not converge, INFO specifies how many
  //                superdiagonals of an intermediate bidiagonal form B
  //                did not converge to zero. See the description of WORK
  //                above for details.
  int INFO = 0;

  // Ready to call the actual factorization routine through PETSc's interface
  LAPACKgesvd_(&JOBU, &JOBVT, &M, &N, &(matrix.get_values()[0]), &LDA, &(sigma_val[0]), &(U_val[0]),
               &LDU, &(VT_val[0]), &LDVT, &(WORK[0]), &LWORK, &INFO);

  // Check return value for errors
  if (INFO != 0)
  {
    Messages::error("during SVD of dense matrix");
  }


  // Load the singular values into sigma, ignore U_val and VT_val
  sigma.resize(sigma_val.size());
  for(unsigned int i=0; i<sigma_val.size(); i++)
    sigma(i) = sigma_val[i];

}
*/
