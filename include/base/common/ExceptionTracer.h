// $Id$

#ifndef _EXCEPTIONTRACER_H_
#define _EXCEPTIONTRACER_H_

#include <execinfo.h>
#include <signal.h>

#include <exception>
#include <iostream>

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
      void* array[25];
      int nSize = backtrace(array, 25);
      char** symbols = backtrace_symbols(array, nSize);

      for (int i = 0; i < nSize; i++)
      {
        std::cerr << symbols[i] << endl;
      }

      free(symbols);
    }
};


#endif // _EXCEPTIONTRACER_H_
