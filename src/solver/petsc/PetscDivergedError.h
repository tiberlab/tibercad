// $Id$

#ifndef _PETSCDIVERGEDERROR_H_
#define _PETSCDIVERGEDERROR_H_

#include "PetscRuntimeError.h"

#include "petsc_macro.h"


/**
 * This exception signals a convergence failure of a PETSc solver
 *
 *
 * Use one of the derived classes instead of this one.
 *
 * The argument @c reason should be one of @c KSPConvergedReason
 * or @c SNESConvergedReason (cf. PETSc documentation).
 */
class PetscDivergedError : public PetscRuntimeError
{

  public:

    enum PETScSolverType { UNKNOWN = 0, KSP, SNES };
    
    PetscDivergedError(int reason, int iteration, double fnorm,
        PETScSolverType type = UNKNOWN);


    int get_iteration(void) const { return _iteration; };
    double get_fnorm(void) const { return _fnorm; };
    int get_reason(void) const { return _reason; };

    PETScSolverType get_solver_type(void) const { return _type; };


  private:

    int _iteration;
    double _fnorm;
    int _reason;
    PETScSolverType _type;
};

#endif // _PETSCDIVERGEDERROR_H_
