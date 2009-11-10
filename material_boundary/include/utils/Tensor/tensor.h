// Copyright (c) 1999-2003, A.H. van den Boogaard

#ifndef TENSOR_H
#define TENSOR_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

using std::sqrt;

// forward declarations
class Tensor2Sym;
class Tensor4DSym;

class Tensor1
{
   private:
      double a[3];

   public:
      Tensor1(){};
      explicit Tensor1(short i){ assert(i==0); a[0]=a[1]=a[2]=0.0; }

      double& operator()(int i)
      {
         assert( i>0 && i<4 );
         return a[i-1];
      }

      const double& operator() (int i) const
      {
         assert( i>0 && i<4 );
         return a[i-1];
      }

      Tensor1& operator+=(const Tensor1& b)
      {
         a[0]+=b.a[0]; a[1]+=b.a[1]; a[2]+=b.a[2];
         return *this;
      }

      Tensor1& operator-=(const Tensor1& b)
      {
         a[0]-=b.a[0]; a[1]-=b.a[1]; a[2]-=b.a[2];
         return *this;
      }

      Tensor1& operator*=(double s)
      {
         a[0]*=s; a[1]*=s; a[2]*=s;
         return *this;
      }

      Tensor1& operator/=(double s)
      {
         double m=1.0/s;
         a[0]*=m; a[1]*=m; a[2]*=m;
         return *this;
      }

      Tensor1& operator=(short i){ assert(i==0); a[0]=a[1]=a[2]=0.0; return *this; }
};


// Virtual base class for second order tensors

class Tensor2
{
   private:

   protected:
      Tensor2(){}

   public:
      virtual ~Tensor2(){}

      virtual const double& operator() (int i, int j) const=0;
      virtual void invariants( double *I, double *II, double *III ) const;
};


// General second order tensor

class Tensor2Gen : public Tensor2
{
   private:
      double a[9];

   public:
      Tensor2Gen(){};
      explicit Tensor2Gen(short i);
      Tensor2Gen( const Tensor2Sym &t );

      Tensor2Gen transpose() const;

      Tensor2Gen& operator=(short i);

      Tensor2Gen& operator+=( const Tensor2Gen &B )
      {
         for ( int i=0; i<9; i++ ) a[i] += B.a[i];
         return *this;
      }

      Tensor2Gen& operator-=(const Tensor2Gen& B)
      {
         for ( int i=0; i<9; i++ ) a[i] -= B.a[i];
         return *this;
      }

      Tensor2Gen& operator*=(double s)
      {
         for ( int i=0; i<9; i++ ) a[i] *= s;
         return *this;
      }

      Tensor2Gen& operator/=(double s)
      {
         double m=1.0/s;
         for ( int i=0; i<9; i++ ) a[i] *= m;
         return *this;
      }

      double& operator()(int i, int j)
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         return a[(j-1)*3+i-1];
      }

      const double& operator() (int i, int j) const
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         return a[(j-1)*3+i-1];
      }

      Tensor1 operator() (int i) const
      {
         // return column i
         Tensor1 r;
         r(1) = (*this)(1,i);
         r(2) = (*this)(2,i);
         r(3) = (*this)(3,i);
         return r;
      }

      friend Tensor2Gen inv( const Tensor2Gen &A );
      friend double doubleContraction( const Tensor2Gen &A, const Tensor2Gen &B );
};


// Symmetric second order tensor

class Tensor2Sym : public Tensor2
{
   // put data only on lower left triangle
 
   private:
      double a[6];  // store by diagonal: 11, 22, 33, 21, 32, 31

      void Givens( int p, int q, Tensor2Gen *V=0 );

   public:
      Tensor2Sym(){};
      explicit Tensor2Sym( short i );

      inline double& operator()(int i, int j)
      {  // only lower left reference (because it is adaptable)
         assert( i>0 && i<4 );
         assert( j>0 && j<=i );
         switch ( i-j ){
            case 0: return a[i-1];
            case 1: return a[i+1];
            case 2: return a[5];
            default: throw(1);
         }
      }

      inline const double& operator() (int i, int j) const
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         if ( j>i ){ int k(j); j=i; i=k; }
         switch ( i-j ){
            case 0: return a[i-1];
            case 1: return a[i+1];
            case 2: return a[5];
            default: throw(1);
         }
      }

      void invariants( double *I, double *II, double *III ) const;
      void eigen( double *lambda1, double *lambda2, double *lambda3,
                  Tensor2Gen *V=0 ) const;

      Tensor2Sym& operator=(short i);

      Tensor2Sym& operator+=( const Tensor2Sym &B )
      {
         for ( int i=0; i<6; i++ ) a[i] += B.a[i];
         return *this;
      }

      Tensor2Sym& operator-=(const Tensor2Sym& B)
      {
         for ( int i=0; i<6; i++ ) a[i] -= B.a[i];
         return *this;
      }

      Tensor2Sym& operator*=(double s)
      {
         for ( int i=0; i<6; i++ ) a[i] *= s;
         return *this;
      }

      Tensor2Sym& operator/=(double s)
      {
         double m=1.0/s;
         for ( int i=0; i<6; i++ ) a[i] *= m;
         return *this;
      }

      void getVector( double *A ) const
         { const double sqrt2( sqrt(2.0) );
           for ( int i=0; i<3; i++ ) A[i] = a[i];
           for ( int i=3; i<6; i++ ) A[i] = sqrt2*a[i]; }
      void putVector( const double *A )
         { const double sqrt2inv( 1.0/sqrt(2.0) );
           for ( int i=0; i<3; i++ ) a[i] = A[i];
           for ( int i=3; i<6; i++ ) a[i] = sqrt2inv*A[i]; }

      friend Tensor2Sym inv( const Tensor2Sym &A );
      friend Tensor2Sym square( const Tensor2Sym &A );
      friend double doubleContraction( const Tensor2Sym &A, const Tensor2Sym &B );
      friend Tensor2Sym doubleContraction( const Tensor4DSym &A, const Tensor2Sym &B );
};


// Virtual base class for fourth order tensors

class Tensor4
{
   private:

   protected:
      Tensor4(){}

   public:
      virtual ~Tensor4(){}

      virtual const double& operator() ( int i, int j, int k, int l ) const=0;
};


// General fourth order tensor

class Tensor4Gen : public Tensor4
{
   private:
      double a[81];

   public:
      Tensor4Gen(){};
      explicit Tensor4Gen( short i );
      Tensor4Gen( const Tensor4DSym &t );

      double& operator()( int i, int j, int k, int l )
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         assert( k>0 && k<4 );
         assert( l>0 && l<4 );
         return a[(i-1)*27+(j-1)*9+(k-1)*3+(l-1)];
      }

      const double& operator() ( int i, int j, int k, int l ) const
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         assert( k>0 && k<4 );
         assert( l>0 && l<4 );
         return a[(i-1)*27+(j-1)*9+(k-1)*3+(l-1)];
      }

      Tensor4Gen& operator=(short i);

      Tensor4Gen& operator+=( const Tensor4Gen &A )
      {
         for ( int i=0; i<81; i++ ) a[i] += A.a[i];
         return *this;
      }

      Tensor4Gen& operator-=( const Tensor4Gen &A )
      {
         for ( int i=0; i<81; i++ ) a[i] -= A.a[i];
         return *this;
      }

      Tensor4Gen& operator*=(double s)
      {
         for ( int i=0; i<81; i++ ) a[i] *= s;
         return *this;
      }

      Tensor4Gen& operator/=(double s)
      {
         double m=1.0/s;
         for ( int i=0; i<81; i++ ) a[i] *= m;
         return *this;
      }
};


// 4th order tensor with symmetry in minor and major index

class Tensor4DSym : public Tensor4
{
   // double symmetric: A_ijkl = A_klij = A_jikl
   // for every symmetry store lower left part (i>j,k>l,ij>kl)

   private:
      double a[21];

   public:
      Tensor4DSym(){};
      explicit Tensor4DSym( short i );

      double& operator()(int i, int j, int k, int l )
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<=i );
         assert( k>0 && k<4 );
         assert( l>0 && l<=k );
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
         assert( I>=J );

         switch ( I-J ){
            case 0: return a[I];
            case 1: return a[I+5];
            case 2: return a[I+9];
            case 3: return a[I+12];
            case 4: return a[I+14];
            case 5: return a[20];
            default: throw(1);
         }
      }

      const double& operator() (int i, int j, int k, int l ) const
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         if ( j>i ){ int ii(j); j=i; i=ii; }
         assert( k>0 && k<4 );
         assert( l>0 && l<4 );
         if ( l>k ){ int kk(l); l=k; k=kk; }
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
         if ( J>I ){ int II(J); J=I; I=II; }

         switch ( I-J ){
            case 0: return a[I];
            case 1: return a[I+5];
            case 2: return a[I+9];
            case 3: return a[I+12];
            case 4: return a[I+14];
            case 5: return a[20];
            default: throw(1);
         }
      }

      Tensor4DSym& operator=(short i);

      Tensor4DSym& operator+=( const Tensor4DSym &A )
      {
         for ( int i=0; i<21; i++ ) a[i] += A.a[i];
         return *this;
      }

      Tensor4DSym& operator-=( const Tensor4DSym &A )
      {
         for ( int i=0; i<21; i++ ) a[i] -= A.a[i];
         return *this;
      }

      Tensor4DSym& operator*=(double s)
      {
         for ( int i=0; i<21; i++ ) a[i] *= s;
         return *this;
      }

      Tensor4DSym& operator/=(double s)
      {
         double m=1.0/s;
         for ( int i=0; i<21; i++ ) a[i] *= m;
         return *this;
      }

      void getFortranMatrix( double *A ) const;
      void putFortranMatrix( const double *A );

      friend Tensor2Sym doubleContraction( const Tensor4DSym &A, const Tensor2Sym &B );
};


// 4th order tensor with symmetry in minor index

class Tensor4LSym : public Tensor4
{
   // lower index symmetric: A_ijkl = A_jikl = A_ijlk
   // for every symmetry store lower left part (i>j,k>l)

   private:
      double a[36];

   public:
      Tensor4LSym(){};
      explicit Tensor4LSym( short i );
      Tensor4LSym( const Tensor4DSym &S );

      double& operator()(int i, int j, int k, int l )
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<=i );
         assert( k>0 && k<4 );
         assert( l>0 && l<=k );
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

         return a[6*J+I];
      }

      const double& operator() (int i, int j, int k, int l ) const
      {
         assert( i>0 && i<4 );
         assert( j>0 && j<4 );
         if ( j>i ){ int ii(j); j=i; i=ii; }
         assert( k>0 && k<4 );
         assert( l>0 && l<4 );
         if ( l>k ){ int kk(l); l=k; k=kk; }
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

         return a[6*J+I];
      }

      Tensor4LSym& operator=(short i);

      Tensor4LSym& operator+=( const Tensor4LSym &A )
      {
         for ( int i=0; i<36; i++ ) a[i] += A.a[i];
         return *this;
      }

      Tensor4LSym& operator-=( const Tensor4LSym &A )
      {
         for ( int i=0; i<36; i++ ) a[i] -= A.a[i];
         return *this;
      }

      Tensor4LSym& operator*=(double s)
      {
         for ( int i=0; i<36; i++ ) a[i] *= s;
         return *this;
      }

      Tensor4LSym& operator/=(double s)
      {
         double m=1.0/s;
         for ( int i=0; i<36; i++ ) a[i] *= m;
         return *this;
      }

      void getFortranMatrix( double *A ) const;
      void putFortranMatrix( const double *A );
};


// output stream operators

std::ostream& operator<<(std::ostream &s, const Tensor1 &T);
std::ostream& operator<<(std::ostream &s, const Tensor2 &T);
std::ostream& operator<<(std::ostream &s, const Tensor4Gen &T);
std::ostream& operator<<(std::ostream &s, const Tensor4DSym &T);
std::ostream& operator<<(std::ostream &s, const Tensor4LSym &T);


// arithmetic operators +, - and scalar * and /

inline Tensor1 operator+(const Tensor1& a, const Tensor1& b)
                        { Tensor1 r(a); return r+=b; }
inline Tensor1 operator-(const Tensor1& a, const Tensor1& b)
                        { Tensor1 r(a); return r-=b; }
inline Tensor1 operator*(const Tensor1& a, double b)
                        { Tensor1 r(a); return r*=b; }
inline Tensor1 operator*(double b, const Tensor1& a)
                        { Tensor1 r(a); return r*=b; }
inline Tensor1 operator/(const Tensor1& a, double b)
                        { Tensor1 r(a); return r/=b; }

inline Tensor2Gen operator+(const Tensor2Gen& A, const Tensor2Gen& B)
                           { Tensor2Gen R(A); return R+=B; }
inline Tensor2Gen operator-(const Tensor2Gen& A, const Tensor2Gen& B)
                           { Tensor2Gen R(A); return R-=B; }
inline Tensor2Gen operator*(const Tensor2Gen& A, double s)
                           { Tensor2Gen R(A); return R*=s; }
inline Tensor2Gen operator*(double s, const Tensor2Gen& A)
                           { Tensor2Gen R(A); return R*=s; }
inline Tensor2Gen operator/(const Tensor2Gen& A, double s)
                           { Tensor2Gen R(A); return R/=s; }

inline Tensor2Sym operator+(const Tensor2Sym& A, const Tensor2Sym& B)
                           { Tensor2Sym R(A); return R+=B; }
inline Tensor2Sym operator-(const Tensor2Sym& A, const Tensor2Sym& B)
                           { Tensor2Sym R(A); return R-=B; }
inline Tensor2Sym operator*(const Tensor2Sym& A, double s)
                           { Tensor2Sym R(A); return R*=s; }
inline Tensor2Sym operator*(double s, const Tensor2Sym& A)
                           { Tensor2Sym R(A); return R*=s; }
inline Tensor2Sym operator/(const Tensor2Sym& A, double s)
                           { Tensor2Sym R(A); return R/=s; }

inline Tensor4Gen operator+(const Tensor4Gen& A, const Tensor4Gen& B)
                           { Tensor4Gen R(A); return R+=B; }
inline Tensor4Gen operator-(const Tensor4Gen& A, const Tensor4Gen& B)
                           { Tensor4Gen R(A); return R-=B; }
inline Tensor4Gen operator*(const Tensor4Gen& A, double s)
                           { Tensor4Gen R(A); return R*=s; }
inline Tensor4Gen operator*(double s, const Tensor4Gen& A)
                           { Tensor4Gen R(A); return R*=s; }
inline Tensor4Gen operator/(const Tensor4Gen& A, double s)
                           { Tensor4Gen R(A); return R/=s; }

inline Tensor4DSym operator+(const Tensor4DSym& A, const Tensor4DSym& B)
                            { Tensor4DSym R(A); return R+=B; }
inline Tensor4DSym operator-(const Tensor4DSym& A, const Tensor4DSym& B)
                            { Tensor4DSym R(A); return R-=B; }
inline Tensor4DSym operator*(const Tensor4DSym& A, double s)
                            { Tensor4DSym R(A); return R*=s; }
inline Tensor4DSym operator*(double s, const Tensor4DSym& A)
                            { Tensor4DSym R(A); return R*=s; }
inline Tensor4DSym operator/(const Tensor4DSym& A, double s)
                            { Tensor4DSym R(A); return R/=s; }

inline Tensor4LSym operator+(const Tensor4LSym& A, const Tensor4LSym& B)
                            { Tensor4LSym R(A); return R+=B; }
inline Tensor4LSym operator-(const Tensor4LSym& A, const Tensor4LSym& B)
                            { Tensor4LSym R(A); return R-=B; }
inline Tensor4LSym operator*(const Tensor4LSym& A, double s)
                            { Tensor4LSym R(A); return R*=s; }
inline Tensor4LSym operator*(double s, const Tensor4LSym& A)
                            { Tensor4LSym R(A); return R*=s; }
inline Tensor4LSym operator/(const Tensor4LSym& A, double s)
                            { Tensor4LSym R(A); return R/=s; }


// unary operators + and -

inline Tensor1 operator-(const Tensor1& a){ return -1.0*a; }
inline Tensor1 operator+(const Tensor1& a){ return a; }
inline Tensor2Gen operator-(const Tensor2Gen& A){ return -1.0*A; }
inline Tensor2Gen operator+(const Tensor2Gen& A){ return A; }
inline Tensor2Sym operator-(const Tensor2Sym& A){ return -1.0*A; }
inline Tensor2Sym operator+(const Tensor2Sym& A){ return A; }
inline Tensor4Gen operator-(const Tensor4Gen& A){ return -1.0*A; }
inline Tensor4Gen operator+(const Tensor4Gen& A){ return A; }
inline Tensor4DSym operator-(const Tensor4DSym& A){ return -1.0*A; }
inline Tensor4DSym operator+(const Tensor4DSym& A){ return A; }
inline Tensor4LSym operator-(const Tensor4LSym& A){ return -1.0*A; }
inline Tensor4LSym operator+(const Tensor4LSym& A){ return A; }


// Tensor-Tensor multiplications

inline double operator*( const Tensor1& a, const Tensor1& b )
{
   return ( a(1)*b(1) + a(2)*b(2) + a(3)*b(3) );
}

inline Tensor1 operator*(const Tensor1& a, const Tensor2& B)
{
   Tensor1 r;
   double a1=a(1), a2=a(2), a3=a(3);

   r(1) = a1*B(1,1)+a2*B(2,1)+a3*B(3,1);
   r(2) = a1*B(1,2)+a2*B(2,2)+a3*B(3,2);
   r(3) = a1*B(1,3)+a2*B(2,3)+a3*B(3,3);

   return r;
}

inline Tensor1 operator*( const Tensor2 &A, const Tensor1 &b )
{
   Tensor1 r;
   double b1=b(1), b2=b(2), b3=b(3);

   r(1) = A(1,1)*b1 + A(1,2)*b2 + A(1,3)*b3;
   r(2) = A(2,1)*b1 + A(2,2)*b2 + A(2,3)*b3;
   r(3) = A(3,1)*b1 + A(3,2)*b2 + A(3,3)*b3;

   return r;
}

inline Tensor2Gen operator*(const Tensor2& A, const Tensor2& B)
{
   Tensor2Gen T2;
   double b;

   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
      {
         b = 0;
         for ( int k=1; k<=3; k++ ) b += A(i,k)*B(k,j);
         T2(i,j) = b;
      }

   return T2;
}


double doubleContraction( const Tensor2 &A, const Tensor2 &B );

Tensor2Gen doubleContraction( const Tensor4 &A, const Tensor2 &B );
Tensor2Gen doubleContraction( const Tensor2 &A, const Tensor4 &B );
Tensor2Sym doubleContraction( const Tensor4LSym &A, const Tensor2Sym &B );

Tensor4Gen doubleContraction( const Tensor4 &A, const Tensor4 &B );


inline Tensor2Gen operator*(const Tensor4Gen& A, const Tensor2Gen& B)
                           {return doubleContraction(A,B); }
inline Tensor2Gen operator*(const Tensor2Gen& A, const Tensor4Gen& B)
                           {return doubleContraction(A,B); }

inline Tensor2Sym operator*(const Tensor4DSym& A, const Tensor2Sym& B)
                           {return doubleContraction(A,B); }
inline Tensor2Sym operator*(const Tensor2Sym& A, const Tensor4DSym& B)
                           {return doubleContraction(B,A); }

inline Tensor2Sym operator*(const Tensor4LSym& A, const Tensor2Sym& B)
                           {return doubleContraction(A,B); }
inline Tensor2Sym operator*(const Tensor2Sym& A, const Tensor4LSym& B)
                           {return doubleContraction(B,A); }


// Specialized multiplications with symmetric result

Tensor2Sym multAtA( const Tensor2 &A ); // A(transpose)*A
Tensor2Sym multAAt( const Tensor2 &A ); // A*A(transpose)
Tensor2Sym multAtSA( const Tensor2 &A, const Tensor2Sym &S ); // A(transpose)*S*A
Tensor2Sym multASAt( const Tensor2 &A, const Tensor2Sym &S ); // A*S*A(transpose)


// Dyadic products

inline Tensor2Gen dyadic( const Tensor1 &a, const Tensor1 &b )
{
   Tensor2Gen T2;
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
         T2(i,j) = a(i)*b(j);
   return T2;
}

inline Tensor2Sym selfDyadic( const Tensor1 &a )
{
   Tensor2Sym T2;
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
         T2(i,j) = a(i)*a(j);
   return T2;
}

inline Tensor4Gen dyadic( const Tensor2 &A, const Tensor2 &B )
{
   Tensor4Gen T4;
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=3; j++ )
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=3; l++ )
               T4( i,j,k,l) = A(i,j)*B(k,l);

   return T4;
}

inline Tensor4LSym dyadic( const Tensor2Sym &A, const Tensor2Sym &B )
{
   Tensor4LSym T4;
   for ( int i=1; i<=3; i++ )
      for ( int j=1; j<=i; j++ )
         for ( int k=1; k<=3; k++ )
            for ( int l=1; l<=k; l++ )
               T4( i,j,k,l) = A(i,j)*B(k,l);

   return T4;
}

inline Tensor4DSym selfDyadic( const Tensor2Sym &A )
{
   Tensor4DSym T4(0);
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

               T4( i,j,k,l) = A(i,j)*A(k,l);
            }
      }

   return T4;
}

Tensor1 vectorProduct( const Tensor1& a, const Tensor1& b );


// Tensor functions

Tensor2Sym sym( const Tensor2 &A );
Tensor2Gen skew( const Tensor2 &A );

Tensor2Gen dev( const Tensor2Gen &A );
Tensor2Sym dev( const Tensor2Sym &A );
Tensor2Sym log( const Tensor2Sym &S );

double norm( const Tensor1 &a );
double norm( const Tensor2 &A );
double trace( const Tensor2 &A );
double det( const Tensor2 &A );

void findroots( const Tensor2Sym &S, Tensor2Sym* U, Tensor2Sym* U_inv );


// Polar decompositions

void polarRight( const Tensor2Gen &F, Tensor2Gen *R, Tensor2Sym *U );
void polarLeft( const Tensor2Gen &F, Tensor2Sym *V, Tensor2Gen *R );


#endif  // TENSOR_H
