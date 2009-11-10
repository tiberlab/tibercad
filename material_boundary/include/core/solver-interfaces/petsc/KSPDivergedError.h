// $Id$

#ifndef _KSPDIVERGEDERROR_H_
#define _KSPDIVERGEDERROR_H_

#include "PetscDivergedError.h"


/**
 * This exception signals a convergence failure of the linear PETSc solver
 *
 * The argument @c reason should be one of @c KSPConvergedReason
 * (cf. PETSc documentation).
 */
class KSPDivergedError : public PetscDivergedError
{

  public:
    KSPDivergedError(int reason, int iteration, double fnorm)
      : PetscDivergedError(reason, iteration, fnorm, KSP) {};

  private:

};

#endif // _KSPDIVERGEDERROR_H_
