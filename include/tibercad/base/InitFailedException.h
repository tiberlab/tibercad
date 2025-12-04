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
 * \file InitFailedException.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef _INITFAILEDEXCEPTION_H_
#define _INITFAILEDEXCEPTION_H_

#include "tibercad/base/ExceptionTracer.h"

#include <stdexcept>
#include <string>

//! An exception class for failed initialisation
class InitFailedException : public std::runtime_error, ExceptionTracer
{

  public:
    InitFailedException(const char* msg)
      : std::runtime_error(msg) {};

    InitFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};

    //void * array[25];
    //size_t entries = backtrace(array, sizeof(array) / sizeof(void*));
    //char ** symbols = backtrace_symbols(array, entries);
    //for ( size_t i = 2; i < entries; i++ ) {
    //   cerr <<  symbols[i] << endl;
    //}


#endif // _INITFAILEDEXCEPTION_H_
