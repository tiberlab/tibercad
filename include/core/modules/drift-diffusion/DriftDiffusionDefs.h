// $Id$

#ifndef _DRIFTDIFFUSIONDEFS_H_
#define _DRIFTDIFFUSIONDEFS_H_


#include "TiberModule.h"


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
    ECURRENT = 0x02,
    HCURRENT = 0x04,
    FULLYCOUPLED = 0x07,
    CURRENTS = 0x06,
    ELECTRONS = 0x02,
    HOLES = 0x04,
    BOTH = 0x06,
    EQUILIBRIUM = 0x10
  };

  //! The known basic Generation/Recombination processes
  enum RecombinationModel
  {
    SRH         = 0x01, //< Shockley-Read-Hall
    AUGER       = 0x02, //< Auger
    DIRECT      = 0x04  //< direct (radiative)
  };

  //! The variables used in Drift-Diffusion
  enum DDVariable
  {
    POTENTIAL = 0, //< electrical potential
    FERMIE    = 1, //< electron electro-chemical potential
    FERMIH    = 2, //< hole electro-chemical potential
    DENSE     = 3, //< electron density
    DENSH     = 4, //< hole density
  };

}


#endif // _DRIFTDIFFUSIONDEFS_H_
