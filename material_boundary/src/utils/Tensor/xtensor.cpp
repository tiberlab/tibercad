// Copyright (c) 2001, A.H. van den Boogaard

#include "xtensor.h"
#include "tutils.h"
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;

Tensor2Sym solve( const Tensor4DSym &A, const Tensor2Sym& rhs )
{
   double D[36];
   double y[6];

   A.getFortranMatrix( D );
   rhs.getVector( y );

   int N(6);

#if 1
   double C[6], Q[6];

   if ( QR_factor( D, &N, C, Q ) !=0 )
   {
      cout << "QR_factor failed."  << endl;
      exit(1);
   }

   if ( QR_solve(D,&N,C,Q,y) )
   {
      cout << "QR_solve did not work." << endl;
      exit(1);
   }
#else
   int piv[6];
   if ( LU_Factor( D, &N, piv )!=0 )
   {
      cout << "LU_Factor failed."  << endl;
      exit(1);
   }
 
   if ( LU_Solve( D, &N, piv, y )!=0 )
   {
      cout << "LU_Solve did not work." << endl;
      exit(1);
   }
#endif

   Tensor2Sym solution;

   solution.putVector( y );

   return solution;
}

Tensor2Sym deviatoric_solve( const Tensor4DSym &A, const Tensor2Sym& rhs )
{
   assert( fabs(rhs(1,1)+rhs(2,2)+rhs(3,3)) < 1e-6*(1+norm(rhs)) );

   double B[36];

   A.getFortranMatrix( B );

   double D[25];
   int k(0);
   for ( int i=0; i<6; i++ )
   {
      if (i==2) continue;
      int l(0);
      for ( int j=0; j<6; j++ )
      {
         if (j==2) continue;
         if (l<2) D[k+5*l] = B[i+6*j]-B[i+6*2];
         else     D[k+5*l] = B[i+6*j];
         l++;
      }
      k++;
   }

   double x[6];
   rhs.getVector( x );

   double y[5];

   y[0] = x[0];
   y[1] = x[1];
   y[2] = x[3];
   y[3] = x[4];
   y[4] = x[5];

   int N(5);

#if 1
   double C[5], Q[5];

   if ( QR_factor( D, &N, C, Q ) !=0 )
   {
      cout << "QR_factor failed."  << endl;
      exit(1);
   }

   if ( QR_solve(D,&N,C,Q,y) )
   {
      cout << "QR_solve did not work." << endl;
      exit(1);
   }
#else
   int piv[5];
   if ( LU_Factor( D, &N, piv )!=0 )
   {
      cout << "LU_Factor failed."  << endl;
      exit(1);
   }
 
   if ( LU_Solve( D, &N, piv, y )!=0 )
   {
      cout << "LU_Solve did not work." << endl;
      exit(1);
   }
#endif

   x[0] = y[0];
   x[1] = y[1];
   x[2] = -y[0]-y[1];
   x[3] = y[2];
   x[4] = y[3];
   x[5] = y[4];

   Tensor2Sym solution;
   solution.putVector( x );

   return solution;
}

Tensor4DSym inv( const Tensor4DSym& A )
{
   double m[36];
   int n=6;
   int p[6];
   double hv[6];

   A.getFortranMatrix( m );

   invFortranMatrix( m, &n, p, hv );

   Tensor4DSym Ainv;
   Ainv.putFortranMatrix( m );

   return( Ainv );
}

Tensor4DSym push_forward( const Tensor4DSym& A, const Tensor2Gen& F )
{
   // calculate B_ijkl = F_ip F_jq F_kr F_ls A_pqrs

   Tensor4DSym B;

   for ( int i=1; i<4; i++ )
      for ( int j=1; j<=i; j++ )
      {
         int I;
         switch ( i-j ){
            case 0: I = i-1;
                    break;
            case 1: I = i+1;
                    break;
            case 2: I = 5;
                    break;
            default:throw(1);
         }

         for ( int k=1; k<4; k++ )
            for ( int l=1; l<=k; l++ )
            {
               int J;
               switch ( k-l ){
                  case 0: J = k-1;
                          break;
                  case 1: J = k+1;
                          break;
                  case 2: J = 5;
                          break;
                  default:throw(1);
               }
               if ( J>I ) continue;

               double value(0);
               for ( int p=1; p<4; p++ )
                  for ( int q=1; q<4; q++ )
                     for ( int r=1; r<4; r++ )
                        for ( int s=1; s<4; s++ )
                           value += F(i,p)*F(j,q)*F(k,r)*F(l,s)*A(p,q,r,s);

               B(i,j,k,l) = value;
            }
      }

   return B;
}


Tensor2Sym solve( const Tensor4LSym &A, const Tensor2Sym& rhs )
{
   double D[36];
   double y[6];

   A.getFortranMatrix( D );
   rhs.getVector( y );

   int N(6);

#if 1
   double C[6], Q[6];

   if ( QR_factor( D, &N, C, Q ) !=0 )
   {
      cout << "QR_factor failed."  << endl;
      exit(1);
   }

   if ( QR_solve(D,&N,C,Q,y) )
   {
      cout << "QR_solve did not work." << endl;
      exit(1);
   }
#else
   int piv[6];
   if ( LU_Factor( D, &N, piv )!=0 )
   {
      cout << "LU_Factor failed."  << endl;
      exit(1);
   }
 
   if ( LU_Solve( D, &N, piv, y )!=0 )
   {
      cout << "LU_Solve did not work." << endl;
      exit(1);
   }
#endif

   Tensor2Sym solution;
   solution.putVector( y );

   return solution;
}

Tensor2Sym deviatoric_solve( const Tensor4LSym &A, const Tensor2Sym& rhs )
{
   assert( fabs(rhs(1,1)+rhs(2,2)+rhs(3,3)) < 1e-6*(1+norm(rhs)) );

   double B[36];
   A.getFortranMatrix( B );

   double D[25];
   int k(0);
   for ( int i=0; i<6; i++ )
   {
      if (i==2) continue;
      int l(0);
      for ( int j=0; j<6; j++ )
      {
         if (j==2) continue;
         if (l<2) D[k+5*l] = B[i+6*j]-B[i+6*2];
         else     D[k+5*l] = B[i+6*j];
         l++;
      }
      k++;
   }

   double x[6];
   rhs.getVector( x );

   double y[5];

   y[0] = x[0];
   y[1] = x[1];
   y[2] = x[3];
   y[3] = x[4];
   y[4] = x[5];

   int N(5);

#if 1
   double C[5], Q[5];

   if ( QR_factor( D, &N, C, Q ) !=0 )
   {
      cout << "QR_factor failed."  << endl;
      exit(1);
   }

   if ( QR_solve(D,&N,C,Q,y) )
   {
      cout << "QR_solve did not work." << endl;
      exit(1);
   }
#else
   int piv[5];
   if ( LU_Factor( D, &N, piv )!=0 )
   {
      cout << "LU_Factor failed."  << endl;
      exit(1);
   }
 
   if ( LU_Solve( D, &N, piv, y )!=0 )
   {
      cout << "LU_Solve did not work." << endl;
      exit(1);
   }
#endif

   Tensor2Sym solution;
   x[0] = y[0];
   x[1] = y[1];
   x[2] = -y[0]-y[1];
   x[3] = y[2];
   x[4] = y[3];
   x[5] = y[4];

   solution.putVector( x );
   return solution;
}

Tensor4LSym inv( const Tensor4LSym& A )
{
   double m[36];
   int n=6;
   int p[6];
   double hv[6];

   A.getFortranMatrix( m );

   invFortranMatrix( m, &n, p, hv );

   Tensor4LSym Ainv;
   Ainv.putFortranMatrix( m );

   return( Ainv );
}

Tensor4LSym push_forward( const Tensor4LSym& A, const Tensor2Gen& F )
{
   // calculate B_ijkl = F_ip F_jq F_kr F_ls A_pqrs

   Tensor4LSym B;

   for ( int i=1; i<4; i++ )
      for ( int j=1; j<=i; j++ )
      {
         for ( int k=1; k<4; k++ )
            for ( int l=1; l<=k; l++ )
            {
               double value(0);
               for ( int p=1; p<4; p++ )
                  for ( int q=1; q<4; q++ )
                     for ( int r=1; r<4; r++ )
                        for ( int s=1; s<4; s++ )
                           value += F(i,p)*F(j,q)*F(k,r)*F(l,s)*A(p,q,r,s);

               B(i,j,k,l) = value;
            }
      }

   return B;
}

Tensor2Gen expW( const Tensor2Gen& W )
{
   // return the exponent of the skew symmetric tensor W
   // usually W is the skew part of the velocity gradient times time increment

   double omega;

   omega = sqrt( 0.5*doubleContraction(W,W) );
   
   Tensor2Gen R(1);
   
   if ( omega > 1e-10 )
   {
      R += sin(omega)/omega * W + (1-cos(omega))/(omega*omega) * W*W;
   }

   return R;
}

