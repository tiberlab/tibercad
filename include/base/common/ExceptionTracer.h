// $Id$

#ifndef _EXCEPTIONTRACER_H_
#define _EXCEPTIONTRACER_H_

#ifdef DEBUG
#include <execinfo.h>
#include <iostream>
#include <cstdlib>
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
      void* array[25];
      int nSize = backtrace(array, 25);
      char** symbols = backtrace_symbols(array, nSize);

      for (int i = 0; i < nSize; i++)
      {
        std::cerr << symbols[i] << std::endl;
      }

      free(symbols);
#endif
    }
};


#endif // _EXCEPTIONTRACER_H_
