// $Id: DriftDiffusionDefs.h 3542 2013-03-01 09:31:59Z maufder $

#ifndef _DRIFTDIFFUSIONDEFS_H_
#define _DRIFTDIFFUSIONDEFS_H_



//! A namespace for drift-diffusion specific definitions
/*!
 */
namespace DriftDiffusionDefs
{

  //! Specifies the type of coupling used
  /*!
   * The values are hardcoded to hex values so they can be easily
   * checked using logical operators
   */
  enum Coupling
  {
    POISSON = 0x01,
    FULLYCOUPLED = 0x07,
    CURRENTS = 0x06,
    EQUILIBRIUM = 0x10
  };

}


#endif // _DRIFTDIFFUSIONDEFS_H_
