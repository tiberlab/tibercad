// $Id$


#ifndef _DEVICEEXCEPTION_H_
#define _DEVICEEXCEPTION_H_

#include "tibercad/base/ExceptionTracer.h"

#include <stdexcept>
#include <string>

//! An exception class for failed initialisation
class DeviceException : public std::runtime_error, ExceptionTracer
{

  public:
    DeviceException(const char* msg)
      : std::runtime_error(msg) {};

    DeviceException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _DEVICEEXCEPTION_H_
