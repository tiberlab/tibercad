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
 * \file SolverException.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_SOLVEREXCEPTION_H
#define TC_SOLVEREXCEPTION_H

#include <stdexcept>
#include <string>

//! An exception class for the solver interfaces
class SolverException : public std::runtime_error
{

  public:
    SolverException(const char* msg)
      : std::runtime_error(msg) {};

    SolverException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // TC_SOLVEREXCEPTION_H
