/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file PetscDivergedError.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_PETSCDIVERGEDERROR_H
#define TC_PETSCDIVERGEDERROR_H

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

#endif // TC_PETSCDIVERGEDERROR_H
