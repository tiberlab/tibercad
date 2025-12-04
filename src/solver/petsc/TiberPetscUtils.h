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
 * \file TiberPetscUtils.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef _TIBERPETSCUTILS_H_
#define _TIBERPETSCUTILS_H_

#include "PetscRuntimeError.h"

#include "petscksp.h"


class ModelOptions;

//! Provides utilities for PETSc
class TiberPetscUtils
{

  public:

    //! Extract KSPType from options
    static const KSPType extract_KSPType(const ModelOptions& options);


    //! Extract PCType from options
    static const PCType extract_PCType(const ModelOptions& options);


    //! Extract line search type (none, cubic, quadratic) from options
    static int extract_LSType(const ModelOptions& options);


    //! Check the PETSc error number
    /*!
     * Throws an exception if errorcode contains a PETSc error number
     */
    static void _checkerr(int errorcode, int line, const char* file);
#ifdef DEBUG
#define checkerr(errcode) _checkerr(errcode, __LINE__, __FILE__);
#else
     static void checkerr(int errorcode)
     {
       if (errorcode != 0)
         throw PetscRuntimeError(errorcode);
     }
#endif


  private:

    TiberPetscUtils();

};


//
// inline members
//



#endif // _TIBERPETSCUTILS_H_
