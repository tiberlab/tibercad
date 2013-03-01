// $Id$

#ifndef _MESSAGES_H_
#define _MESSAGES_H_

#include "tiber_dll.h"

#include <string>
#include <iostream>
#include <sstream>
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
    static void info(const std::string& msg, bool newline = true);


    //! Add empty line
    static void newline(void);

    //! Add a break line (e.g. "<<<<<<<<...")
    static void frameline(const std::string& start, const char c,
        const std::string& name = "");

    //! Set the log file
    static void set_log_file(const std::string& logfile);

    //! Close the log file
    static void close_log_file(void);


    //! Our own endline
    static const std::string endl;


    //! The available text width
    static int available_width(void);


    //! Print statistics on errors and warnings
    static void print_statistics(void);


    //! Tell if we are interactive
    static bool& interactive(void);


    //! Tell if we are should stop at warnings
    static bool& stop_on_warning(void);


  private:


    //! The log file
    static std::ofstream _log;

    //! The maximum line width
    static const int _max_width;

    //! The global indentation level
    static int _indent;

    //! The indentation width
    static const int _indent_width;

    //! The number of errors
    static int _error_count;

    //! The number of warnings
    static int _warning_count;

    //! Interactive or not
    /*!
     * If interactive errors need to be
     * acknowledged by the user.
     */
    static bool _interactive;

    //! Stop also at warnings
    /*!
     * If set, warnings have to be acknowledged by the user
     */
    static bool _stop_on_warning;


    //! The local indentation
    int _indent_loc;
};



//
// inline methods
//

inline
Messages::Messages(void) : _indent_loc(0) { }


inline
Messages::~Messages(void)
{
  //newline();
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


inline
bool&
Messages::interactive(void)
{
  return _interactive;
}

inline
bool&
Messages::stop_on_warning(void)
{
  return _stop_on_warning;
}


#endif // _MESSAGES_H_
