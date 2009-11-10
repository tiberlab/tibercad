// $Id$


#ifndef _SOLVEREXCEPTION_H_
#define _SOLVEREXCEPTION_H_

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



#endif // _SOLVEREXCEPTION_H_
