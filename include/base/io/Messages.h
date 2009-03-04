// $Id$

#ifndef _MESSAGES_H_
#define _MESSAGES_H_

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

    //! Constructor
    Messages(void);

    //! Destructor
    /*!
     * Resets indentation to former value
     */
    ~Messages(void);


    //! Add a level of indentation
    void indent(void);

    //! Decrease level of indentation
    void unindent(void);


    //! Print a warning
    static void warning(const std::string& msg);

    
    //! Print a debug message
    static void debug(const std::string& msg);

    
    //! Print an error
    static void error(const std::string& msg);

    
    //! Print an info
    static void info(const std::string& msg);


    //! Add empty line
    static void newline(void);

    
  private:


    //! Warning keyword
    static std::string _warning;


    //! Error keyword
    static std::string _error;


    //! To reset output format
    static std::string _reset;


    //! An 'extended' endl
    static std::string _endl;


    //! The global indentation level
    static int _indent;

    //! The local indentation
    int _indent_loc;

};



//
// inline methods
// 

inline
Messages::Messages(void) : _indent_loc(0) { };


inline
Messages::~Messages(void)
{
  newline();
  _indent -= _indent_loc;
}


inline
void
Messages::indent(void)
{
  _indent_loc++;
  _indent++;
}


inline
void
Messages::unindent(void)
{
  if (_indent_loc > 0)
  {
    _indent_loc--;
    _indent--;
  }
}



#endif // _MESSAGES_H_
