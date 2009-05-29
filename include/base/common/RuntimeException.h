// $Id$


#ifndef _RUNTIMEEXCEPTION_H_
#define _RUNTIMEEXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for generic runtime errors
class RuntimeException : public std::runtime_error
{

  public:
    RuntimeException(const char* msg)
      : std::runtime_error(msg) {};

    RuntimeException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _RUNTIMEEXCEPTION_H_
