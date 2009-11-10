// $Id$


#ifndef _MODELERROREXCEPTION_H_
#define _MODELERROREXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for failed initialisation
class ModelErrorException : public std::runtime_error
{

  public:
    ModelErrorException(const char* msg)
      : std::runtime_error(msg) {};

    ModelErrorException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _MODELERROREXCEPTION_H_
