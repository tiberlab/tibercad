// $Id$

#ifndef _TIBERCAD_H_
#define _TIBERCAD_H_

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
  extern char** cmdline_argv;

  //! We remember the number of command line arguments
  extern int cmdline_argc;


  //! An initialization routine
  /*!
   * This routine calls init() of libmesh and other libraries, if needed
   */
  void init(int argc, char** argv);


  //! A cleanup routine
  /*!
   * This routine calls close() of libmesh and other libraries, if needed
   */
  int cleanup(void);

}


#endif // _TIBERCAD_H_
