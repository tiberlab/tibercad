#ifndef _OPTGPL_H_
#define _OPTGPL_H_

#include "Vff.h"

extern "C"
{
#include "cg_user.h"
#include "cg_descent.h"
}
/*! This is a static library used to call optimization methods
 *  The code is based on a C GPL library :
 *
 *    |      A conjugate gradient method with guaranteed descent       |
      |             C-code Version 1.1  (October 6, 2005)              |
      |                    Version 1.2  (November 14, 2005)            |
      |                    Version 2.0  (September 23, 2007)           |
      |                    Version 3.0  (May 18, 2008)                 |
      |                    Version 4.0  (March 28, 2011)               |
      |                    Version 4.1  (April 8, 2011)                |
      |                    Version 4.2  (April 14, 2011)               |
      |                    Version 5.0  (May 1, 2011)                  |
      |                                                                |
      |           William W. Hager    and   Hongchao Zhang             |
      |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |                 Copyright by William W. Hager                  |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |________________________________________________________________|
 *
 *  and has been adapted to be compatible with ANSI C++ compiler
 *  IMPORTANT: the library IS NOT thread safe
 */
class OptGpl
  {
  public:

  OptGpl(Vff& vff);

  ~OptGpl(void);

  void solve(double tolerance = 1e-5);

  private:

  Vff& _vff;

   static int _n_istances;

   int cg_descent(
       double            *x, /* input: starting guess, output: the solution */
       INT                n, /* problem dimension */
       cg_stats       *Stat, /* structure with statistics (can be NULL) */
       cg_parameter  *UParm, /* user parameters, NULL = use default parameters */
       double      grad_tol, /* StopRule = 1: |g|_infty <= max (grad_tol,
                                              StopFac*initial |g|_infty) [default]
                                StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
       double         *Work  /* either size 4n work array or NULL */
       );

   int cg_evaluate
   (
     char    *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
     char     *nan, /* y means check function/derivative values for nan */
     cg_com   *Com
   );

   double cg_cubic
   (
     double  a,
     double fa, /* function value at a */
     double da, /* derivative at a */
     double  b,
     double fb, /* function value at b */
     double db  /* derivative at b */
   );

   void cg_copy
   (
     double *y, /* output of copy */
     double *x, /* input of copy */
     int     n  /* length of vectors */
   );

   void cg_step
   (
     double *xtemp, /*output vector */
     double     *x, /* initial vector */
     double     *d, /* search direction */
     double  alpha, /* stepsize */
     INT         n  /* length of the vectors */
   );

   void cg_default
   (
     cg_parameter   *Parm
   );

   void cg_printParms
   (
     cg_parameter  *Parm
   );

   double cg_dot
   (
     double *x, /* first vector */
     double *y, /* second vector */
     INT     n /* length of vectors */
   );

   int cg_contract
   (
     double    *A, /* left side of bracketing interval */
     double   *fA, /* function value at a */
     double   *dA, /* derivative at a */
     double    *B, /* right side of bracketing interval */
     double   *fB, /* function value at b */
     double   *dB, /* derivative at b */
     cg_com  *Com  /* cg com structure */
   );

   int cg_line
   (
     cg_com   *Com /* cg com structure */
   );

   int cg_tol
   (
     double     gnorm, /* gradient sup-norm */
     cg_com    *Com    /* cg com */
   );

   int cg_Wolfe
   (
     double   alpha, /* stepsize */
     double       f, /* function value associated with stepsize alpha */
     double    dphi, /* derivative value associated with stepsize alpha */
     cg_com    *Com  /* cg com */
   );


  };


#endif
