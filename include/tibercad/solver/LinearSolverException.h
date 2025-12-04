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
 * \file LinearSolverException.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef _LINEARSOLVEREXCEPTION_H_
#define _LINEARSOLVEREXCEPTION_H_

#include "tibercad/solver/SolverException.h"

#include <stdexcept>
#include <string>

//! An exception class for the linear solver interfaces
class LinearSolverException : public SolverException
{

  public:
    LinearSolverException(const char* msg)
      : SolverException(msg) {};

    LinearSolverException(const std::string& msg)
      : SolverException(msg) {};


  private:

};



#endif // _LINEARSOLVEREXCEPTION_H_
