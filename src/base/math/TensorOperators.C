// $Id$

#include "TensorOperators.h"



RealTensorValue doubleContraction(const Tensor4DSym& A, const RealTensorValue& B)
{
  /*
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
   */

  RealTensorValue T(B);
  RealTensorValue R;

  return R;
}
