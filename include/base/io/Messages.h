// $Id$

#ifndef _MESSAGES_H_
#define _MESSAGES_H_

#include "TypeDefs.h"

#include <string>
#include <iostream>

#ifdef error
# undef error
#endif

//! Print messages to standard out or standard error
/*!
 * This class contains methods to print messages to stderr and stdout
 * depending on compilation mode (debug, optimized, profiling)
 */
class Messages
{

  public:

    //! Print a warning
    static void warning(const std::string& msg);

    
    //! Print a debug message
    static void debug(const std::string& msg);

    
    //! Print an error
    static void error(const std::string& msg);

    
    //! Print an info
    static void info(const std::string& msg);


    
  private:

    //! This class is for static use only!
    Messages(void);


    //! Warning keyword
    static std::string _warning;


    //! Error keyword
    static std::string _error;


    //! To reset output format
    static std::string _reset;


    //! An 'extended' endl
    static std::string _endl;

};



//
// inline methods
// 


inline
void
Messages::debug(const std::string& msg)
{
#ifdef DEBUG
  error(msg);
#else
  ignore_unused_variable(msg);
#endif
}




#endif // _MESSAGES_H_
