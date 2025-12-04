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
 * \file PetscRuntimeError.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _PETSCRUNTIMEERROR_H_
#define _PETSCRUNTIMEERROR_H_

#include "tibercad/solver/SolverException.h"

#include <string>

class PetscRuntimeError : public SolverException
{

  public:

    PetscRuntimeError(int reason);

    virtual ~PetscRuntimeError(void) throw() {};

    int get_reason(void) const { return _reason; };

    virtual const char* what(void) const throw();

  protected:

    void set_message(const std::string& msg);

  private:

    std::string _msg;

    int _reason;
};


#endif // _PETSCRUNTIMEERROR_H_
