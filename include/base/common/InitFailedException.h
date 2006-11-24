// $Id$


#ifndef _INITFAILEDEXCEPTION_H_
#define _INITFAILEDEXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for failed initialisation
class InitFailedException : public std::runtime_error
{

  public:
    InitFailedException(const char* msg)
      : std::runtime_error(msg) {};

    InitFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _INITFAILEDEXCEPTION_H_
