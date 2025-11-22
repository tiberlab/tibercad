// $Id$


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
