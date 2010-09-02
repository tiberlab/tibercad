// $Id$

#ifndef _TIBERCAD_H_
#define _TIBERCAD_H_

#include "tiber_dll.h"

#include <string>

class Control;


//! Some useful common definitions for TiberCAD
namespace TiberCad 
{

  //! The statistics to be used
  enum Statistics
  {
    BOLTZMANN,     /*!< Boltzmann statistics */
    FERMIDIRAC,    /*!< Fermi-Dirac statistics */
    BOSEEINSTEIN   /*!< Bose-Einstein statistics */
  };


  //! Types of symmetries
  enum Symmetry
  {
    NONE,         /*!< no special symmetry */
    CYLINDRICAL   /*!< cylinder symmetry */
  };


  //! We keep a pointer to the command line arguments
  extern char** cmdline_argv TBDLLOCAL;

  //! We remember the number of command line arguments
  extern int cmdline_argc TBDLLOCAL;


  //! The installation root directory
  /*!
   * This is read from the environment
   */
  extern std::string tiberroot;


  //! An initialization routine
  /*!
   * This routine calls init() of libmesh and other libraries, if needed.
   */
  void init(int argc, char** argv);


  //! A cleanup routine
  /*!
   * This routine calls close() of libmesh and other libraries, if needed
   */
  int cleanup(void);


  //! The TiberCAD major version
  extern const int MajorVersion;


  //! The TiberCAD minor version
  extern const int MinorVersion;


  //! The TiberCAD subminor version
  extern const int SubMinorVersion;


  //! The TiberCAD subversion release number
  extern const int SvnRevision;


  //! The full TiberCAD version string
  std::string TiberCadVersion(bool include_svn_release = true);


  //! Get access to the control module
  Control& get_control(void);


  //! Get the directory where to put output files
  const std::string& get_output_dir(void);


  //! Get the complete suffix for the output filenames
  const std::string& get_filename_suffix(void);


  //! Clear the suffix for the output filenames
  void clear_filename_suffix(void);


  //! Append something to the suffix for the output filenames
  /*!
   * The filename suffix will be appended to all output files which
   * contain plot data.
   * The suffix itself will be prepended by a '_'
   */
  void append_to_filename_suffix(const std::string& suffix);


  //! Prepend something to the suffix for the output filenames
  /*!
   * The filename suffix will be appended to all output files which
   * contain plot data.
   * The suffix itself will be prepended by a '_'
   */
  void prepend_to_filename_suffix(const std::string& suffix);


  //! Delete the first output filename suffix part
  void drop_first_filename_suffix(void);


  //! Delete the last output filename suffix part
  void drop_last_filename_suffix(void);



  //! Get the output format
  /*!
   * \return a string that identifies the type of output files
   * to generate
   *
   * Currently the following formats are supported:
   * \li \c gmv for GMV
   * \li \c ise for Tecplot
   * \li \c gnu for GnuPlot
   */
  const std::string& get_output_format(void);

}


#endif // _TIBERCAD_H_
