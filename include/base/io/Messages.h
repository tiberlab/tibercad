// $Id$

#ifndef _MESSAGES_H_
#define _MESSAGES_H_

#include <string>
#include <iostream>
#include <fstream>

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


    //! Set the log file
    static void set_log_file(const std::string& logfile);

    //! Close the log file
    static void close_log_file(void);


    //! Our own endline
    static const std::string endl;


    //! The available text width
    static int available_width(void);

  private:


    //! The log file
    static std::ofstream _log;

    //! The maximum line width
    static const int _max_width;

    //! The global indentation level
    static int _indent;

    //! The indentation width
    static const int _indent_width;

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


inline
int
Messages::available_width(void)
{
  return _max_width - _indent * _indent_width;
}


#endif // _MESSAGES_H_
