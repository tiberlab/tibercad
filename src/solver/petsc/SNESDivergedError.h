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
 * \file SNESDivergedError.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_SNESDIVERGEDERROR_H
#define TC_SNESDIVERGEDERROR_H

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

#endif // TC_SNESDIVERGEDERROR_H
