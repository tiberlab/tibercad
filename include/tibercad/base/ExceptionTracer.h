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
 * \file ExceptionTracer.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _EXCEPTIONTRACER_H_
#define _EXCEPTIONTRACER_H_


#ifdef DEBUG 
#if !defined(__CYGWIN__) && !defined(__MINGW32__)
#include <execinfo.h>
#include <iostream>
#include <cstdlib>
#endif
#endif

//! A class for exception tracing
/*!
 * This code is inspired by the article
 * "C++ exception-handling tricks for Linux" by Sachin Agrawal
 * (http://www.ibm.com/developerworks/linux/library/l-cppexcep.html)
 */
class ExceptionTracer
{

  public:

    ExceptionTracer(void)
    {
#ifdef DEBUG
#if !defined(__CYGWIN__) && !defined(__MINGW32__)
      void* array[25];
      int nSize = backtrace(array, 25);
      char** symbols = backtrace_symbols(array, nSize);

      for (int i = 0; i < nSize; i++)
      {
        std::cerr << symbols[i] << std::endl;
      }

      free(symbols);
#endif
#endif
    }
};


#endif // _EXCEPTIONTRACER_H_
