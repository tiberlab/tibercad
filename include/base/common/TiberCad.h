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

}


#endif // _TIBERCAD_H_
