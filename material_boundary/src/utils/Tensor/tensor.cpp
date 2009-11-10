// Copyright (c) 1999-2003, A.H. van den Boogaard

#include "tensor.h"
#include <iomanip>
#include <cmath>

#ifdef NO_STD_LIMITS
#include <float.h>
#else
#include <limits>
#endif

// unnamed namespace: local to this file
namespace
{
#ifdef NO_STD_LIMITS
const double ZERO = 1000*DBL_MIN;
#else
const double ZERO = 1000*std::numeric_limits<double>::min();
#endif

const double ZERO_TOL=1.e-12;
}

using std::ostream;
using std::setw;
using std::endl;

// Tensor2 member functions

void Tensor2::invariants( double *I, double *II, double *III ) const
{
   if ( I  !=0 ) *I   = trace( *this );
   if ( II !=0 ) *II  = (*this)(1,1)*(*this)(2,2)-(*this)(2,1)*(*this)(1,2) +
                        (*this)(1,1)*(*this)(3,3)-(*this)(3,1)*(*this)(1,3) +
                        (*this)(2,2)*(*this)(3,3)-(*this)(3,2)*(*this)(2,3);
   if ( III!=0 ) *III = det( *this );
}


// Tensor2Gen constructors and assignment

Tensor2Gen::Tensor2Gen( const Tensor2Sym &t )
{
   a[0] = t(1,1);
   a[1] = a[3] = t(2,1);
   a[2] = a[6] = t(3,1);
   a[4] = t(2,2);
   a[5] = a[7] = t(3,2);
   a[8] = t(3,3);
}

Tensor2Gen::Tensor2Gen(short index)
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   for ( int j=0; j<9; j++ ) a[j]=0;
   if ( index==1 ) a[0] = a[4] = a[8] = 1;
}

Tensor2Gen& Tensor2Gen::operator=(short index)
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   for ( int j=0; j<9; j++ ) a[j]=0;
   if ( index==1 ) a[0] = a[4] = a[8] = 1;
   return *this;
}


// return transpose of this tensor

Tensor2Gen Tensor2Gen::transpose() const
{
   Tensor2Gen R;
   R.a[0] = a[0]; R.a[1] = a[3]; R.a[2] = a[6];
   R.a[3] = a[1]; R.a[4] = a[4]; R.a[5] = a[7];
   R.a[6] = a[2]; R.a[7] = a[5]; R.a[8] = a[8];
   return R;
}


// constructor and assigment of Tensor2Sym

Tensor2Sym::Tensor2Sym(short index)
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<6; j++ ) a[j]=0;
             break;
     case 1: for ( int j=0; j<3; j++ ) a[j]=1;
             for ( int j=3; j<6; j++ ) a[j]=0;
   }
}

Tensor2Sym& Tensor2Sym::operator=(short index)
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<6; j++ ) a[j]=0;
             break;
     case 1: for ( int j=0; j<3; j++ ) a[j]=1;
             for ( int j=3; j<6; j++ ) a[j]=0;
   }
   return *this;
}


// Perform Givens rotation of this symmetric tensor

void Tensor2Sym::Givens( int p, int q, Tensor2Gen *V )
{ 
   /* Givens transformation for row and column p and q
      See Matrix Computations, Golub & Van Loan, 2nd ed.,
      section 5.1.8 and section 8.5.2 */

   double Apq, tau, t, s, c, taup, tauq;
   Tensor2Sym &A(*this);

   assert( p>0 && p<4 );
   assert( q>0 && q<4 );
   assert( p!=q );

   if ( p<q ) // swap p and q
   {
      int r(p);
      p = q;
      q = r;
   }

   Apq = A(p,q);
   if ( fabs(Apq) > ZERO )
   {
      double App = a[p-1]; //A(p,p)
      double Aqq = a[q-1]; //A(q,q)

      tau = (Aqq - App)/(2*Apq);
      t = 1/(std::fabs(tau)+sqrt(1+tau*tau));
      if ( tau < 0 ) t = -t;
      c = 1/sqrt(1+t*t);
      s = t*c;

      double Apq2cs = Apq*2*c*s;
      double c2 = c*c;
      double s2 = s*s;

      // pre multiply by J.transpose() and post multiply by J

      a[p-1] = App*c2 + Aqq*s2 - Apq2cs; // A(p,p) 
      a[q-1] = Aqq*c2 + App*s2 + Apq2cs; // A(q,q) 
      A(p,q) = 0;

      if ( q==1 )
         if ( p==2 )
         {
            double Arp = a[4]; //A(3,2)
            double Arq = a[5]; //A(3,1)
            a[4] = Arp*c-Arq*s;
            a[5] = Arp*s+Arq*c;
         }
         else // p==3
         {
            double Arp = a[4]; //A(3,2)
            double Arq = a[3]; //A(2,1)
            a[4] = Arp*c-Arq*s;
            a[3] = Arp*s+Arq*c;
         }
      else // p==3 && q==2
      {
            double Arp = a[5]; //A(3,1)
            double Arq = a[3]; //A(2,1)
            a[5] = Arp*c-Arq*s;
            a[3] = Arp*s+Arq*c;
      }

      // post multiply V by J to update eigen vectors (as columns)
      if ( V != 0 )
      {
         taup = (*V)(1,p); tauq = (*V)(1,q);
         (*V)(1,p) = c*taup-s*tauq; (*V)(1,q) = s*taup+c*tauq;
         taup = (*V)(2,p); tauq = (*V)(2,q);
         (*V)(2,p) = c*taup-s*tauq; (*V)(2,q) = s*taup+c*tauq;
         taup = (*V)(3,p); tauq = (*V)(3,q);
         (*V)(3,p) = c*taup-s*tauq; (*V)(3,q) = s*taup+c*tauq;
      }
   }
}


// Determine the eigenvalues and optionally eigenvectors of this symmetric
// tensor. The eigenvectors are the columns of tensor V

void Tensor2Sym::eigen( double *lambda1, double *lambda2, double *lambda3,
                      Tensor2Gen *V ) const
{
   Tensor2Sym T2(*this);

   if ( V!=0 ) *V=1;

   if ( T2.a[3]*T2.a[3]+T2.a[4]*T2.a[4]+T2.a[5]*T2.a[5] > ZERO )
   {
      double delta = 1.e-20*doubleContraction(T2,T2);//pow(1.e-10*norm(T2),2);
      int iter = 0;
      do
      {
         /* For large diagonal terms and small but nonzero off-diagonals,
            do at least one sweep */
         T2.Givens(2,1,V);
         T2.Givens(3,1,V);
         T2.Givens(3,2,V);
      } while ( T2.a[3]*T2.a[3]+T2.a[4]*T2.a[4]+T2.a[5]*T2.a[5] > delta
             && ++iter < 5 );
   }

   *lambda1 = T2.a[0];
   *lambda2 = T2.a[1];
   *lambda3 = T2.a[2];
}

void Tensor2Sym::invariants( double *I, double *II, double *III ) const
{
   if ( I  !=0 ) *I   = a[0]+a[1]+a[2];
   if ( II !=0 ) *II  = a[0]*a[1] + a[0]*a[2] + a[1]*a[2]
                       -a[3]*a[3] - a[4]*a[4] - a[5]*a[5];
   if ( III!=0 ) *III = a[0]*(a[1]*a[2]-a[4]*a[4]) - a[3]*(a[3]*a[2]-a[4]*a[5])
                       +a[5]*(a[3]*a[4]-a[1]*a[5]);
}


// Constructors and assignment of a general 4th order tensor

Tensor4Gen::Tensor4Gen( const Tensor4DSym &S )
{
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
               (*this)(i,j,k,l) = S(i,j,k,l);
}

Tensor4Gen::Tensor4Gen( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<81; j++ ) a[j]=0;
             break;
     case 1: for ( int i=1; i<=3; i++ )
                for ( int j=1; j<=3; j++ )
                   for ( int k=1; k<=3; k++ )
                      for ( int l=1; l<=3; l++ )
                         (*this)(i,j,k,l) = (i==k && j==l)?1:0;
   }
}

Tensor4Gen& Tensor4Gen::operator=( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<81; j++ ) a[j]=0;
             break;
     case 1: for ( int i=1; i<=3; i++ )
                for ( int j=1; j<=3; j++ )
                   for ( int k=1; k<=3; k++ )
                      for ( int l=1; l<=3; l++ )
                         (*this)(i,j,k,l) = (i==k && j==l)?1:0;
   }
   return *this;
}


// Constructor and assignment of a symmetric 4th order tensor

Tensor4DSym::Tensor4DSym( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<21; j++ ) a[j]=0;
             break;
     case 1: // 0.5*( delta_ik delta_jl + delta_jk delta_il )
             a[0] = a[1] = a[2] = 1.0;
             a[3] = a[4] = a[5] = 0.5;
             for ( int j=6; j<21; j++ ) a[j] = 0.0;
   }
}

Tensor4DSym& Tensor4DSym::operator=( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<21; j++ ) a[j]=0;
             break;
     case 1: // 0.5*( delta_ik delta_jl + delta_jk delta_il )
             a[0] = a[1] = a[2] = 1.0;
             a[3] = a[4] = a[5] = 0.5;
             for ( int j=6; j<21; j++ ) a[j] = 0.0;
   }
   return *this;
}

// Fill a columnwise stored double array e.g. to use linear algebra libraries

void Tensor4DSym::getFortranMatrix( double *A ) const
{
   const double sqrt2( sqrt(2.0) );

   for ( int i=0; i<6; i++ )
   {
      for ( int j=0; j<=i; j++ )
      {
         double value;
         switch ( i-j ){
            case 0: value = a[i];
                    break;
            case 1: value = a[i+5];
                    break;
            case 2: value = a[i+9];
                    break;
            case 3: value = a[i+12];
                    break;
            case 4: value = a[i+14];
                    break;
            case 5: value = a[20];
                    break;
            default:throw(1);
         }
         if ( i>2 )
         {
            if ( j<3 ) value *= sqrt2;
            else       value *= 2;
         }
         A[i+6*j] = value;
         if ( i!=j ) A[j+6*i] = value;
      }
   }
}

// Get data from a columnwise stored double array
// only the lower-left part is used

void Tensor4DSym::putFortranMatrix( const double *A )
{
   const double sqrt2inv( 1.0/sqrt(2.0) );

   for ( int i=0; i<6; i++ )
   {
      for ( int j=0; j<=i; j++ )
      {
         int k;
         switch ( i-j ){
            case 0: k = i;
                    break;
            case 1: k = i+5;
                    break;
            case 2: k = i+9;
                    break;
            case 3: k = i+12;
                    break;
            case 4: k = i+14;
                    break;
            case 5: k = 20;
                    break;
            default:throw(1);
         }
         if ( i<3 )
            a[k] = A[i+6*j];
         else if ( j<3 )
            a[k] = sqrt2inv*A[i+6*j];
         else
            a[k] = 0.5*A[i+6*j];
      }
   }
}


// Constructor and assignment of a lower symmetric 4th order tensor

Tensor4LSym::Tensor4LSym( const Tensor4DSym &S )
{
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=k; l++ )
               (*this)(i,j,k,l) = S(i,j,k,l);
}

Tensor4LSym::Tensor4LSym( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<36; j++ ) a[j]=0;
             break;
     case 1: // 0.5*( delta_ik delta_jl + delta_jk delta_il )
             for ( int j=1; j<35; j++ ) a[j] = 0.0;
             a[0] = a[7] = a[14] = 1.0;
             a[21] = a[28] = a[35] = 0.5;
   }
}

Tensor4LSym& Tensor4LSym::operator=( short index )
{
   assert( index==0 || index==1 );  // 0: zero tensor, 1: unit tensor

   switch ( index ){
     case 0: for ( int j=0; j<21; j++ ) a[j]=0;
             break;
     case 1: // 0.5*( delta_ik delta_jl + delta_jk delta_il )
             for ( int j=1; j<35; j++ ) a[j] = 0.0;
             a[0] = a[7] = a[14] = 1.0;
             a[21] = a[28] = a[35] = 0.5;
   }
   return *this;
}

// Fill a columnwise stored double array e.g. to use linear algebra libraries

void Tensor4LSym::getFortranMatrix( double *A ) const
{
   for ( int i=0; i<36; i++ ) A[i] = a[i];

   const double sqrt2( sqrt(2.0) );

   for ( int i=0; i<3; i++ )
      for ( int j=3; j<6; j++ )
      {
         A[i+6*j] *= sqrt2;
         A[j+6*i] *= sqrt2;
      }
   for ( int i=3; i<6; i++ )
      for ( int j=3; j<6; j++ )
         A[i+6*j] *= 2.0;
}

// Get data from a columnwise stored double array

void Tensor4LSym::putFortranMatrix( const double *A )
{
   const double sqrt2inv( 1.0/sqrt(2.0) );

   int k(0);
   for ( int i=0; i<6; i++ )
      for ( int j=0; j<6; j++, k++ )
      {
         if ( i<3 )
         {
            if ( j<3 )
               a[k] = A[k];
            else
               a[k] = A[k]*sqrt2inv;
         }
         else
         {
            if ( j<3 )
               a[k] = A[k]*sqrt2inv;
            else
               a[k] = A[k]*0.5;
         }
      }
}



// Non-member functions


// output stream operators

ostream& operator<<(ostream &s, const Tensor1 &T)
{
  int w = s.width();
  s << setw(w) << T(1) << setw(w) << T(2) << setw(w) << T(3) << endl;
  return s;
}

ostream& operator<<(ostream &s, const Tensor2 &T)
{
  int w = s.width();
  s << setw(w) << T(1,1) << setw(w) << T(1,2) << setw(w) << T(1,3) << endl;
  s << setw(w) << T(2,1) << setw(w) << T(2,2) << setw(w) << T(2,3) << endl;
  s << setw(w) << T(3,1) << setw(w) << T(3,2) << setw(w) << T(3,3) << endl;
  return s;
}

ostream& operator<<(ostream &s, const Tensor4Gen &T)
{
  int w = s.width();

  s << setw(w) << T(1,1,1,1) << setw(w) << T(1,1,1,2) << setw(w) << T(1,1,1,3)
    << setw(w) << T(1,1,2,1) << setw(w) << T(1,1,2,2) << setw(w) << T(1,1,2,3)
    << setw(w) << T(1,1,3,1) << setw(w) << T(1,1,3,2) << setw(w) << T(1,1,3,3) << endl;
  s << setw(w) << T(1,2,1,1) << setw(w) << T(1,2,1,2) << setw(w) << T(1,2,1,3)
    << setw(w) << T(1,2,2,1) << setw(w) << T(1,2,2,2) << setw(w) << T(1,2,2,3)
    << setw(w) << T(1,2,3,1) << setw(w) << T(1,2,3,2) << setw(w) << T(1,2,3,3) << endl;
  s << setw(w) << T(1,3,1,1) << setw(w) << T(1,3,1,2) << setw(w) << T(1,3,1,3)
    << setw(w) << T(1,3,2,1) << setw(w) << T(1,3,2,2) << setw(w) << T(1,3,2,3)
    << setw(w) << T(1,3,3,1) << setw(w) << T(1,3,3,2) << setw(w) << T(1,3,3,3) << endl;
  s << setw(w) << T(2,1,1,1) << setw(w) << T(2,1,1,2) << setw(w) << T(2,1,1,3)
    << setw(w) << T(2,1,2,1) << setw(w) << T(2,1,2,2) << setw(w) << T(2,1,2,3)
    << setw(w) << T(2,1,3,1) << setw(w) << T(2,1,3,2) << setw(w) << T(2,1,3,3) << endl;
  s << setw(w) << T(2,2,1,1) << setw(w) << T(2,2,1,2) << setw(w) << T(2,2,1,3)
    << setw(w) << T(2,2,2,1) << setw(w) << T(2,2,2,2) << setw(w) << T(2,2,2,3)
    << setw(w) << T(2,2,3,1) << setw(w) << T(2,2,3,2) << setw(w) << T(2,2,3,3) << endl;
  s << setw(w) << T(2,3,1,1) << setw(w) << T(2,3,1,2) << setw(w) << T(2,3,1,3)
    << setw(w) << T(2,3,2,1) << setw(w) << T(2,3,2,2) << setw(w) << T(2,3,2,3)
    << setw(w) << T(2,3,3,1) << setw(w) << T(2,3,3,2) << setw(w) << T(2,3,3,3) << endl;
  s << setw(w) << T(3,1,1,1) << setw(w) << T(3,1,1,2) << setw(w) << T(3,1,1,3)
    << setw(w) << T(3,1,2,1) << setw(w) << T(3,1,2,2) << setw(w) << T(3,1,2,3)
    << setw(w) << T(3,1,3,1) << setw(w) << T(3,1,3,2) << setw(w) << T(3,1,3,3) << endl;
  s << setw(w) << T(3,2,1,1) << setw(w) << T(3,2,1,2) << setw(w) << T(3,2,1,3)
    << setw(w) << T(3,2,2,1) << setw(w) << T(3,2,2,2) << setw(w) << T(3,2,2,3)
    << setw(w) << T(3,2,3,1) << setw(w) << T(3,2,3,2) << setw(w) << T(3,2,3,3) << endl;
  s << setw(w) << T(3,3,1,1) << setw(w) << T(3,3,1,2) << setw(w) << T(3,3,1,3)
    << setw(w) << T(3,3,2,1) << setw(w) << T(3,3,2,2) << setw(w) << T(3,3,2,3)
    << setw(w) << T(3,3,3,1) << setw(w) << T(3,3,3,2) << setw(w) << T(3,3,3,3) << endl;
  return s;
}

ostream& operator<<(ostream &s, const Tensor4LSym &T)
{
  int w = s.width();

  s << setw(w) << T(1,1,1,1) << setw(w) << T(1,1,2,2) << setw(w) << T(1,1,3,3)
    << setw(w) << T(1,1,1,2) << setw(w) << T(1,1,2,3) << setw(w) << T(1,1,3,1) << endl;
  s << setw(w) << T(2,2,1,1) << setw(w) << T(2,2,2,2) << setw(w) << T(2,2,3,3)
    << setw(w) << T(2,2,1,2) << setw(w) << T(2,2,2,3) << setw(w) << T(2,2,3,1) << endl;
  s << setw(w) << T(3,3,1,1) << setw(w) << T(3,3,2,2) << setw(w) << T(3,3,3,3)
    << setw(w) << T(3,3,1,2) << setw(w) << T(3,3,2,3) << setw(w) << T(3,3,3,1) << endl;
  s << setw(w) << T(1,2,1,1) << setw(w) << T(1,2,2,2) << setw(w) << T(1,2,3,3)
    << setw(w) << T(1,2,1,2) << setw(w) << T(1,2,2,3) << setw(w) << T(1,2,3,1) << endl;
  s << setw(w) << T(2,3,1,1) << setw(w) << T(2,3,2,2) << setw(w) << T(2,3,3,3)
    << setw(w) << T(2,3,1,2) << setw(w) << T(2,3,2,3) << setw(w) << T(2,3,3,1) << endl;
  s << setw(w) << T(3,1,1,1) << setw(w) << T(3,1,2,2) << setw(w) << T(3,1,3,3)
    << setw(w) << T(3,1,1,2) << setw(w) << T(3,1,2,3) << setw(w) << T(3,1,3,1) << endl;
  return s;
}

ostream& operator<<(ostream &s, const Tensor4DSym &T)
{
  int w = s.width();

  s << setw(w) << T(1,1,1,1) << setw(w) << T(1,1,2,2) << setw(w) << T(1,1,3,3)
    << setw(w) << T(1,1,1,2) << setw(w) << T(1,1,2,3) << setw(w) << T(1,1,3,1) << endl;
  s << setw(w) << T(2,2,1,1) << setw(w) << T(2,2,2,2) << setw(w) << T(2,2,3,3)
    << setw(w) << T(2,2,1,2) << setw(w) << T(2,2,2,3) << setw(w) << T(2,2,3,1) << endl;
  s << setw(w) << T(3,3,1,1) << setw(w) << T(3,3,2,2) << setw(w) << T(3,3,3,3)
    << setw(w) << T(3,3,1,2) << setw(w) << T(3,3,2,3) << setw(w) << T(3,3,3,1) << endl;
  s << setw(w) << T(1,2,1,1) << setw(w) << T(1,2,2,2) << setw(w) << T(1,2,3,3)
    << setw(w) << T(1,2,1,2) << setw(w) << T(1,2,2,3) << setw(w) << T(1,2,3,1) << endl;
  s << setw(w) << T(2,3,1,1) << setw(w) << T(2,3,2,2) << setw(w) << T(2,3,3,3)
    << setw(w) << T(2,3,1,2) << setw(w) << T(2,3,2,3) << setw(w) << T(2,3,3,1) << endl;
  s << setw(w) << T(3,1,1,1) << setw(w) << T(3,1,2,2) << setw(w) << T(3,1,3,3)
    << setw(w) << T(3,1,1,2) << setw(w) << T(3,1,2,3) << setw(w) << T(3,1,3,1) << endl;
  return s;
}


// Determine the root and the inverse root of a symmetric tensor.
// Algorithm of L.P. Franca, Computers Math. Applic. Vol18 pp 459-466. 1989
// Implementation by Tibor Fulop, 2003.

void findroots( const Tensor2Sym& S, Tensor2Sym *U, Tensor2Sym *U_inv )
{
   double I, II, III;
   S.invariants(&I,&II,&III);
   double k   =  I*I - 3.*II;
 
   if ( k < ZERO_TOL )
   { // isotropic
      double lambda = std::sqrt( I/3 );
      (*U)     = lambda*Tensor2Sym(1);
      (*U_inv) = Tensor2Sym(1)/lambda;
   }
   else
   {
      // The largest eigenvalue
      double l = std::pow(I,3) - 4.5 *I*II + 13.5*III;

      double fi = std::acos( l / std::pow(k, 1.5) );

      double lambda_sqr = ( I+2.*std::sqrt(k)*std::cos(fi/3.) ) / 3.;
      double lambda     = std::sqrt( lambda_sqr );
    
      // Invariants of U
      double III_U = std::sqrt( III );
      double I_U   = lambda + std::sqrt( -lambda_sqr + I + 2*III_U / lambda );
      double II_U  = 0.5*(I_U*I_U - I);
  
      // (*U) = (I_U*III_U*Tensor2Sym(1)+(I_U*I_U-II_U)*S-square(S))
      //        /(I_U*II_U-III_U);

      *U = 1;
      *U *= I_U*III_U;
      *U += (I_U*I_U-II_U)*S;
      *U -= square(S);
      *U /= I_U*II_U-III_U;

      if ( U_inv != 0 )
      {
         // (*U_inv) = (II_U*Tensor2Sym(1) - I_U*(*U) + S ) / III_U;

         *U_inv = 1;
         *U_inv *= II_U;
         *U_inv -= I_U*(*U);
         *U_inv += S;
         *U_inv /= III_U;
      }
   }
}

// Polar decompositions

void polarLeft( const Tensor2Gen &F, Tensor2Sym *V, Tensor2Gen *R )
{
   /* VR decomposition of a regular 2nd order tensor */

   Tensor2Sym B, Vinv;

   B = multAAt( F ); // = V*V
   findroots(B,V,&Vinv);

   *R = Vinv*F;
}

void polarRight( const Tensor2Gen &F, Tensor2Gen *R, Tensor2Sym *U )
{
   /* RU decomposition of a regular 2nd order tensor */

   Tensor2Sym C, Uinv;

   C = multAtA( F );
   findroots(C,U,&Uinv);

   *R = F*Uinv;
}


// Tensor values functions of 1 tensor

Tensor2Sym sym( const Tensor2 &A )
{
   Tensor2Sym R;
   for ( int i=1; i<=3; i++ )
      for( int j=1; j<=i; j++ )
         if ( i==j )
            R(i,j) = A(i,j);
         else
            R(i,j) = 0.5*(A(i,j)+A(j,i));

   return R;
}

Tensor2Gen skew( const Tensor2 &A )
{
   Tensor2Gen R;
   for ( int i=1; i<=3; i++ )
      for( int j=1; j<=3; j++ )
         if ( i==j )
            R(i,j) = 0;
         else
            R(i,j) = 0.5*(A(i,j)-A(j,i));

   return R;
}

Tensor2Gen dev( const Tensor2Gen &A )
{
   Tensor2Gen R(A);
   double p(-trace(A)/3);
   R(1,1) += p;
   R(2,2) += p;
   R(3,3) += p;

   return R;
}

Tensor2Gen inv( const Tensor2Gen &A )
{
   double Sub11, Sub12, Sub13;
   Sub11 = A.a[4]*A.a[8] - A.a[5]*A.a[7];
   Sub12 = A.a[3]*A.a[8] - A.a[6]*A.a[5];
   Sub13 = A.a[3]*A.a[7] - A.a[6]*A.a[4];
   double detA = A.a[0]*Sub11 - A.a[1]*Sub12 + A.a[2]*Sub13;
   assert( std::fabs(detA)>1.e-15 );

   Tensor2Gen Ainv;

   Ainv.a[0] =  Sub11/detA;
   Ainv.a[1] = (-A.a[1]*A.a[8]+A.a[2]*A.a[7])/detA;
   Ainv.a[2] = ( A.a[1]*A.a[5]-A.a[2]*A.a[4])/detA;
   Ainv.a[3] = -Sub12/detA;
   Ainv.a[4] = ( A.a[0]*A.a[8]-A.a[2]*A.a[6])/detA;
   Ainv.a[5] = (-A.a[0]*A.a[5]+A.a[2]*A.a[3])/detA;
   Ainv.a[6] =  Sub13/detA;
   Ainv.a[7] = (-A.a[0]*A.a[7]+A.a[1]*A.a[6])/detA;
   Ainv.a[8] = ( A.a[0]*A.a[4]-A.a[1]*A.a[3])/detA;

   return Ainv;
}

Tensor2Sym dev( const Tensor2Sym &A )
{
   Tensor2Sym R(A);
   double p(-trace(A)/3);
   R(1,1) += p;
   R(2,2) += p;
   R(3,3) += p;

   return R;
}

Tensor2Sym inv( const Tensor2Sym &A )
{
   double Sub11, Sub12, Sub13;
   Sub11 = A.a[1]*A.a[2] - A.a[4]*A.a[4];
   Sub12 = A.a[3]*A.a[2] - A.a[5]*A.a[4];
   Sub13 = A.a[3]*A.a[4] - A.a[5]*A.a[1];
   double detA = A.a[0]*Sub11 - A.a[3]*Sub12 + A.a[5]*Sub13;
   assert( std::fabs(detA)>1.e-15 );

   Tensor2Sym Ainv;

   Ainv.a[0] =  Sub11/detA;
   Ainv.a[1] = ( A.a[0]*A.a[2]-A.a[5]*A.a[5])/detA;
   Ainv.a[2] = ( A.a[0]*A.a[1]-A.a[3]*A.a[3])/detA;
   Ainv.a[3] = -Sub12/detA;
   Ainv.a[4] = (-A.a[0]*A.a[4]+A.a[5]*A.a[3])/detA;
   Ainv.a[5] =  Sub13/detA;

   return Ainv;
}

Tensor2Sym square( const Tensor2Sym &A )
{
   const double *a(A.a);
   double sqr12(a[3]*a[3]);
   double sqr23(a[4]*a[4]);
   double sqr13(a[5]*a[5]);

   Tensor2Sym Asqr;

   Asqr.a[0] = a[0]*a[0] + sqr12 + sqr13;
   Asqr.a[1] = sqr12 + a[1]*a[1] + sqr23;
   Asqr.a[2] = sqr13 + sqr23 + a[2]*a[2];
   Asqr.a[3] = a[3]*(a[0]+a[1])+a[4]*a[5];
   Asqr.a[4] = a[4]*(a[1]+a[2])+a[3]*a[5];
   Asqr.a[5] = a[5]*(a[0]+a[2])+a[3]*a[4];

   return Asqr;
}

double norm( const Tensor1 &a )
{
    double sqn(0);
    for ( size_t i=1; i<4; i++ )
        sqn += a(i)*a(i);

    return sqrt(sqn);
}

double norm( const Tensor2 &A )
{
    double sqn(0);
    for ( size_t i=1; i<4; i++ )
       for ( size_t j=1; j<4; j++ )
          sqn += A(i,j)*A(i,j);

    return sqrt(sqn);
}

double trace( const Tensor2 &A )
{
   return A(1,1)+A(2,2)+A(3,3);
}

double det( const Tensor2 &A )
{
   return  A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
          -A(2,1)*(A(1,2)*A(3,3)-A(3,2)*A(1,3))
          +A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3));
}


Tensor2Sym multAtA( const Tensor2 &A )
{
   Tensor2Sym S;
   double b;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
      {
         b = 0;
         for ( int k=1; k<=3; k++ ) b += A(k,i)*A(k,j);
         S(i,j) = b;
      }

   return S;
}

Tensor2Sym multAAt( const Tensor2 &A )
{
   Tensor2Sym S;
   double b;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
      {
         b = 0;
         for ( int k=1; k<=3; k++ ) b += A(i,k)*A(j,k);
         S(i,j) = b;
      }

   return S;
}

Tensor2Sym multAtSA( const Tensor2 &A, const Tensor2Sym &S )
{
   // A(transpose)*S*A

   Tensor2Gen SA = S*A;
   Tensor2Sym R;
   double b;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
      {
         b = 0;
         for ( int k=1; k<=3; k++ ) b += A(k,i)*SA(k,j);
         R(i,j) = b;
      }

   return R;
}

Tensor2Sym multASAt( const Tensor2 &A, const Tensor2Sym &S )
{
   // A*S*A(transpose)

   Tensor2Gen AS = A*S;
   Tensor2Sym R;
   double b;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
      {
         b = 0;
         for ( int k=1; k<=3; k++ ) b += AS(i,k)*A(j,k);
         R(i,j) = b;
      }

   return R;
}

Tensor2Sym log( const Tensor2Sym& S )
{
   /* determine logarithm of the symmetric tensor S */

   Tensor1 v1, v2, v3;
   Tensor2Gen V;
   double lambda1, lambda2, lambda3;

   S.eigen( &lambda1, &lambda2, &lambda3, &V );

   lambda1 = S(1,1);
   lambda2 = S(2,2);
   lambda3 = S(3,3);

   v1(1) = V(1,1); v1(2) = V(2,1); v1(3) = V(3,1);
   v2(1) = V(1,2); v2(2) = V(2,2); v2(3) = V(3,2);
   v3(1) = V(1,3); v3(2) = V(2,3); v3(3) = V(3,3);

   // reuse of V
   V = std::log(lambda1)*dyadic(v1,v1) + 
       std::log(lambda2)*dyadic(v2,v2) + 
       std::log(lambda3)*dyadic(v3,v3);

   return sym(V);
}


// double contractions and vector products

double doubleContraction( const Tensor2 &A, const Tensor2 &B )
{
   double d=0;

   for ( int i=1; i<4; i++ )
      for ( int j=1; j<4; j++ )
         d += A(i,j)*B(i,j);

   return d;
}

double doubleContraction( const Tensor2Gen &A, const Tensor2Gen &B )
{
   double d=0;

   for ( int i=0; i<9; i++ )
      d += A.a[i]*B.a[i];

   return d;
}

double doubleContraction( const Tensor2Sym &A, const Tensor2Sym &B )
{
   double d=0;

   for ( int i=0; i<3; i++ )
      d += A.a[i]*B.a[i];

   for ( int i=3; i<6; i++ )
      d += 2*A.a[i]*B.a[i];

   return d;
}

Tensor2Gen doubleContraction( const Tensor4 &A, const Tensor2 &B )
{
   Tensor2Gen R;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
      {
         double r=0;
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
               r += A(i,j,k,l)*B(k,l);
         R(i,j) = r;
      }
   return R;
}

Tensor2Gen doubleContraction( const Tensor2 &A, const Tensor4 &B )
{
   Tensor2Gen R;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
      {
         double r=0;
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
               r += A(k,l)*B(k,l,i,j);
         R(i,j) = r;
      }
   return R;
}

Tensor2Sym doubleContraction( const Tensor4DSym &A, const Tensor2Sym &B )
{
   // with these symmetries A*B = B*A

   Tensor2Sym BB(B);
   BB.a[3] *= 2;
   BB.a[4] *= 2;
   BB.a[5] *= 2;

   Tensor2Sym R;

   const double *Aj = A.a;
   const double *bjp=BB.a, *bjm;
   double *rjp=R.a, *rjm;

   for ( int i=0; i<6; i++ )
      *rjp++ = (*Aj++)*(*bjp++);

   for ( int i=1; i<6; i++ )
   {
      rjp=&R.a[i];
      rjm=&R.a[0];
      bjp=&BB.a[i];
      bjm=&BB.a[0];

      for ( int j=i; j<6; j++ )
      {
         *rjm++ += *Aj*(*bjp++);
         *rjp++ += *Aj*(*bjm++);
         Aj++;
      }
   }

   return R;
}

Tensor2Sym doubleContraction( const Tensor4LSym &A, const Tensor2Sym &B )
{
   // with these symmetries A*B != B*A

   Tensor2Sym R;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
      {
         double r=0;
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
               r += A(i,j,k,l)*B(k,l);
         R(i,j) = r;
      }

   return R;
}

Tensor4Gen doubleContraction( const Tensor4 &A, const Tensor4 &B )
{
   // R_ijkl = A_ijmn B_mnkl

   Tensor4Gen R;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
            {
               double r=0;
               for ( int m=1; m<=3; m++ )
                  for ( int n=1; n<=3; n++ )
                     r += A(i,j,m,n)*B(m,n,k,l);
               R(i,j,k,l) = r;
            }
   return R;
}

Tensor1 vectorProduct( const Tensor1& a, const Tensor1& b )
{
   Tensor1 r;

   r(1) =  a(2)*b(3)-a(3)*b(2);
   r(2) = -a(1)*b(3)+a(3)*b(1);
   r(3) =  a(1)*b(2)-a(2)*b(1);

   return r;
}
