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
 * \file DatabaseException.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef _DATABASEEXCEPTION_H_
#define _DATABASEEXCEPTION_H_

#include "tibercad/base/ExceptionTracer.h"

#include <stdexcept>
#include <string>

//! An exception class for failed Database operations
class DatabaseException : public std::runtime_error, ExceptionTracer
{

  public:
    DatabaseException(const char* msg)
      : std::runtime_error(msg) {};

    DatabaseException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _DATABASEEXCEPTION_H_
