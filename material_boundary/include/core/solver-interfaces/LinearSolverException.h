// $Id$


#ifndef _LINEARSOLVEREXCEPTION_H_
#define _LINEARSOLVEREXCEPTION_H_

#include "SolverException.h"

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
