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
 * \file PardisoSolverException.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef TC_PARDISOSOLVEREXCEPTION_H
#define TC_PARDISOSOLVEREXCEPTION_H

#include "tibercad/solver/SolverException.h"

#include <stdexcept>
#include <string>

//! An exception class for the PARDISO solver
class PardisoSolverException : public SolverException
{

 public:
  
  PardisoSolverException(int error)
    : SolverException("Error in Pardiso solver") { };
    


  private:

};





#endif // TC_PARDISOSOLVEREXCEPTION_H
