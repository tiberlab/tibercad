#include "tutils.h"
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

//
// The QR and LU decomposition and related solution of a linear system is
// based on the template algorithm from TNT.
//
// Copied and adapted from file: lu.h and qr.h
// which is part of:
//    Template Numerical Toolkit (TNT) for Linear Algebra
//    see http://math.nist.gov/tnt for updates
//
// "un-templated", specialized for basic C arrays instead of TNT matrix
// and vector objects.
// The following comments are the original comments, based on TNT objects.
//
// Classical QR factorization example, based on Stewart[1973].
//
// This algorithm computes the factorization of a matrix A
// into a product of an orthognal matrix (Q) and an upper triangular 
// matrix (R), such that QR = A.
//
// Parameters:
//
//  A   (in):       Matrix(1:N, 1:N)
//
//  Q   (output):   Matrix(1:N, 1:N), collection of Householder
//                      column vectors Q1, Q2, ... QN
//
//  R   (output):   upper triangular Matrix(1:N, 1:N)
//
// Returns:  
//
//  0 if successful, 1 if A is detected to be singular
//

#include <cmath>      //for sqrt() & fabs()

// Classical QR factorization, based on Stewart[1973].
//
//
// This algorithm computes the factorization of a matrix A
// into a product of an orthognal matrix (Q) and an upper triangular 
// matrix (R), such that QR = A.
//
// Parameters:
//
//  A   (in/out):  On input, A is square, Matrix(1:N, 1:N), that represents
//                  the matrix to be factored.
//
//                 On output, Q and R is encoded in the same Matrix(1:N,1:N)
//                 in the following manner:
//
//                  R is contained in the upper triangular section of A,
//                  except that R's main diagonal is in D.  The lower
//                  triangular section of A represents Q, where each
//                  column j is the vector  Qj = I - uj*uj'/pi_j.
//
//  C  (output):    vector of Pi[j]
//  D  (output):    main diagonal of R, i.e. D(i) is R(i,i)
//
// Returns:  
//
//  0 if successful, 1 if A is detected to be singular
//

int QR_factor(double *A, const int *N, double *C, double *D)
{
   int M(*N);

    int i,j,k;
    double eta, sigma, sum;

    for (k=0; k<M-1; k++)
    {
        // eta = max |M(i,k)|,  for k <= i <= n
        //
        eta = 0;
        for (i=k; i<M; i++)
        {
            double absA = fabs(A[i+M*k]);
            eta = ( absA >  eta ? absA : eta ); 
        }

        if (eta == 0)           // matrix is singular
            return 1;

        // form Qk and premultiply M by it
        //
        for(i=k; i<M; i++)
            A[i+M*k]  = A[i+M*k] / eta;

        sum = 0;
        for (i=k; i<M; i++)
            sum = sum + A[i+M*k]*A[i+M*k];
//        sigma = sign(A[k+M*k]) *  sqrt(sum);
        sigma = (A[k+M*k]>0) ? sqrt(sum) : -sqrt(sum);


        A[k+M*k] = A[k+M*k] + sigma;
        C[k] = sigma * A[k+M*k];
        D[k] = -eta * sigma;

        for (j=k+1; j<M; j++)
        {
            sum = 0;
            for (i=k; i<M; i++)
                sum = sum + A[i+M*k]*A[i+M*j];
            sum = sum / C[k];

            for (i=k; i<M; i++)
                A[i+M*j] = A[i+M*j] - sum * A[i+M*k];
        }

        D[M-1] = A[M*M-1];
        C[M-1] = 0;
    }

    return 0;
}

// modified form of upper triangular solve, except that the main diagonal
// of R (upper portion of A) is in D.
//
int R_solve( const double *A, const int *M, const double *D, double *b )
{
    int i,j;
    int N = *M;

    double sum;

    if (D[N-1] == 0)
        return 1;

    b[N-1] = b[N-1] / D[N-1];

    for (i=N-2; i>=0; i--)
    {
        if (D[i] == 0)
            return 1;
        sum = 0;
        for (j=i+1; j<N; j++)
            sum = sum + A[i+N*j]*b[j];
        b[i] = ( b[i] - sum ) / D[i];
    }

    return 0;
}


int QR_solve(const double *A, const int *M, const double *c,
             const double *d, double *b)
{
    int N=*M;

    int i,j;
    double sum, tau;

    for (j=0; j<N-1; j++)
    {
       // form Q'*b
       sum = 0;
       for (i=j; i<N; i++)
           sum = sum + A[i+N*j]*b[i];
       if (c[j] == 0)
           return 1;
       tau = sum / c[j];
       for (i=j; i<N; i++)
          b[i] = b[i] - tau * A[i+N*j];
    }
    return R_solve(A, M, d, b);        // solve Rx = Q'b
}


// Solve system of linear equations Ax = b.
//
//  Typical usage:
//
//    Matrix(double) A;
//    Vector(Subscript) ipiv;
//    Vector(double) b;
//
//    1)  LU_Factor(A,ipiv);
//    2)  LU_Solve(A,ipiv,b);
//
//   Now b has the solution x.  Note that both A and b
//   are overwritten.  If these values need to be preserved, 
//   one can make temporary copies, as in 
//
//    O)  Matrix(double) T = A;
//    1)  LU_Factor(T,ipiv);
//    1a) Vector(double) x=b;
//    2)  LU_Solve(T,ipiv,x);
//
//   See details below.
//


// right-looking LU factorization algorithm (unblocked)
//
//   Factors matrix A into lower and upper  triangular matrices 
//   (L and U respectively) in solving the linear equation Ax=b.  
//
//
// Args:
//
// A        (input/output) Matrix(1:n, 1:n)  In input, matrix to be
//                  factored.  On output, overwritten with lower and 
//                  upper triangular factors.
//
// indx     (output) Vector(1:n)    Pivot vector. Describes how
//                  the rows of A were reordered to increase
//                  numerical stability.
//
// Return value:
//
// int      (0 if successful, 1 otherwise)
//
//


int LU_Factor( double *A, const int *M, int *indx )
{
    int N=*M;

    if ( N==0 ) return 0;

    int i,j,k;
    int jp;

    double t;

    for (j=0; j<N; j++)
    {

        // find pivot in column j and  test for singularity.

        jp = j;
        t = fabs(A[j+N*j]);
        for (i=j+1; i<N; i++)
            if ( fabs(A[i+N*j]) > t)
            {
                jp = i;
                t = fabs(A[i+N*j]);
            }

        indx[j] = jp;

        // jp now has the index of maximum element 
        // of column j, below the diagonal

        if ( A[jp+N*j] == 0 )                 
            return 1;       // factorization failed because of zero pivot


        if (jp != j)            // swap rows j and jp
            for (k=0; k<N; k++)
            {
                t = A[j+N*k];
                A[j+N*k] = A[jp+N*k];
                A[jp+N*k] =t;
            }

        if (j<N-1)              // compute elements j+1:N of jth column
        {
            // note A(j,j), was A(jp,p) previously which was
            // guaranteed not to be zero (Label #1)
            //
            double recp = 1.0 / A[j+N*j];

            for (k=j+1; k<N; k++)
                A[k+N*j] *= recp;
        }


        if (j<N)
        {
            // rank-1 update to trailing submatrix:   E = E - x*y;
            //
            // E is the region A(j+1:N, j+1:N)
            // x is the column vector A(j+1:N,j)
            // y is row vector A(j,j+1:N)

            int ii,jj;

            for (ii=j+1; ii<N; ii++)
                for (jj=j+1; jj<N; jj++)
                    A[ii+N*jj] -= A[ii+N*j]*A[j+N*jj];
        }
    }

    return 0;
}   


int LU_Solve( const double *A, const int *M, const int *indx, double *b )
{
    int N=*M;
    int    i,ii=-1,ip,j;
    double sum;

    for (i=0;i<N;i++) 
    {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii!=-1)
            for (j=ii;j<=i-1;j++) 
                sum -= A[i+N*j]*b[j];
        else if (sum) ii=i;
            b[i]=sum;
    }

    for (i=N-1;i>=0;i--) 
    {
        sum=b[i];
        for (j=i+1;j<N;j++) 
            sum -= A[i+N*j]*b[j];
        b[i]=sum/A[i+N*i];
    }

    return 0;
}


// driver functions (not from TNT)

void solveQR( double *A, const int *N, double *C, double *D, double *b )
{
   // A[N*N] input matrix (columnwise), output decomposed matrix
   // C[N]   scratch
   // D[N]   scratch
   // b[N]   input right-hand-side, output solution

   if ( QR_factor( A, N, C, D ) )
   {
      cout << "QR_factor failed."  << endl;
      exit(1);
   }

   if ( QR_solve( A, N, C, D, b ) )
   {
      cout << "QR_solve did not work." << endl;
      exit(1);
   }
}

void solveLU( double *A, const int *N, int *piv, double *b )
{
   // A[N*N] input matrix (columnwise), output decomposed matrix
   // p[N]   pivots (scratch)
   // b[N]   input right-hand-side, output solution

   if ( LU_Factor( A, N, piv ) )
   {
      cout << "LU_Factor failed."  << endl;
      exit(1);
   }
 
   LU_Solve( A, N, piv, b );
}

void invFortranMatrix( double *m, int *N, int *p, double *hv )
{
// invert a square matrix of size N in fortran format
//    m(N,N)
//    p(N)
//    hv(N)  

// See page 171 of 'Introduction to numerical analysis'
// by J. Stoer & R. Bulirsch

   int i,j,k,r;

   int n=*N;
   double max,tmp, hlp;
   int tmpint;

   for (i=0; i<n; i++) p[i]=i;

   for (j=0; j<n; j++)
   {
      max=fabs(m[j+n*j]);
      r=j;
      for (i=j+1; i<n; i++) {
         if ((fabs(tmp=m[i+n*j])) > max) {
            max=tmp;
            r=i;
         }
      }

      if (max==0) throw "determinant zero in calculation of inverse matrix";

      if (r>j) {
         for (k=0; k<n; k++) {
            tmp = m[j+n*k];
            m[j+n*k] = m[r+n*k];
            m[r+n*k] = tmp;
         }
         tmpint = p[j];
         p[j] = p[r];
         p[r] = tmpint;
      }

      tmp = 1/m[j+n*j];
      for (i=0; i<n; i++) {
         m[i+n*j] = tmp * m[i+n*j];
      }
      m[j+n*j] = tmp;

      for (k=0; k<n; k++) {
         if (k==j) continue;
         for (i=0; i<n; i++) {
            if (i==j) continue;
            hlp = m[i+n*k] - (m[i+n*j] * m[j+n*k]);
            m[i+n*k] = hlp;
         }
         m[j+n*k] = -tmp * m[j+n*k];
      }
   }

   for (i=0; i<n; i++) {
      for (k=0; k<n; k++) {
         hv[p[k]]=m[i+n*k];
      }
      for (k=0; k<n; k++) {
         m[i+n*k] = hv[k];
      }
   }
}
