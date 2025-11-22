// $Id$

#ifndef _SNESDIVERGEDERROR_H_
#define _SNESDIVERGEDERROR_H_

#include "PetscDivergedError.h"


/**
 * This exception signals a convergence failure of the nonlinear
 * PETSc solver
 *
 * The argument @c reason should be one of @c SNESConvergedReason
 * (cf. PETSc documentation).
 */
class SNESDivergedError : public PetscDivergedError
{

  public:
    SNESDivergedError(int reason, int iteration, double fnorm)
      : PetscDivergedError(reason, iteration, fnorm, SNES) {};

  private:

};

#endif // _SNESDIVERGEDERROR_H_
