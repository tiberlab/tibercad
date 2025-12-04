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
 * \file PetscRuntimeError.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "PetscRuntimeError.h"

#include <sstream>


PetscRuntimeError::PetscRuntimeError(int reason)
  : SolverException("Internal PETSc error."),
    _reason(reason)
{
  std::ostringstream os;
  os << "Internal PETSc error: " << get_reason();
  _msg = os.str();
}


const char*
PetscRuntimeError::what(void) const throw()
{
  return _msg.c_str();
}


void
PetscRuntimeError::set_message(const std::string& msg)
{
  _msg = msg;
}
