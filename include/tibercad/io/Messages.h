/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file Messages.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_MESSAGES_H
#define TC_MESSAGES_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#ifdef error
# undef error
#endif

namespace libMesh {
  namespace Parallel {
    class Communicator;
  }
}


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


    //! Print an hint
    static void hint(const std::string& msg, bool newline = true);


    //! Add empty line
    static void newline(void);

    //! Add a break line (e.g. "<<<<<<<<...")
    static void frameline(const std::string& start, const char c,
        const std::string& name = "");

    //! Set the log file
    static void set_log_file(const std::string& logfile,
        libMesh::Parallel::Communicator& comm, int rank);

    //! Close the log file
    static void close_log_file(void);

    //! Set the communicator
    static void set_communicator(
        libMesh::Parallel::Communicator& comm, int rank);

    //! Set the stdout stream
    static void set_stdout(std::ostream& os = nullstream);


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

    //! The stdout
    static std::ostream* _cout;

    //! A stream writing to /dev/null
    static std::ofstream nullstream;

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


    //! The MPI communicator associated with the logger
    static libMesh::Parallel::Communicator* _mpi_comm;

    //! The rank on the communicator that should write to the file
    static int _rank;
};



//
// inline methods
//




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


#endif // TC_MESSAGES_H
