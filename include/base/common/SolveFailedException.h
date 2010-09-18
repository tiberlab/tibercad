// $Id$


#ifndef _SOLVEFAILEDEXCEPTION_H_
#define _SOLVEFAILEDEXCEPTION_H_

#include "ExceptionTracer.h"

#include <stdexcept>
#include <string>

//! An exception class for failed solve
class SolveFailedException : public std::runtime_error, ExceptionTracer
{

  public:
    SolveFailedException(const char* msg)
      : std::runtime_error(msg) {};

    SolveFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _SOLVEFAILEDEXCEPTION_H_
