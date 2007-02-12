// $Id$

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

/**
 * Physical constants
 *
 * SI units, converted to cm as length unit
 *
 */
namespace Constants
{
  /**
   * Boltzmann constant in eV / K
   */
  extern const double k_Boltzmann;
  extern const double& k_B;
  extern const double& kb;

  /**
   * Elementary charge in C
   */
  extern const double elementary_charge;
  extern const double& e;

  /**
   * 1 eV
   */
  extern const double electron_volt;
  extern const double& eV;

  /**
   * dielectric constant of vacuum
   */
  extern const double epsilon;
  extern const double& eps;
  extern const double& e0;

  /**
   * rest mass of an electron in kg
   */
  extern const double electron_mass;
  extern const double& me;

  /**
   * Planck's constant in eV * s
   */
  extern const double plancks_constant;
  extern const double& h;
  extern const double hbar;


  
  //!Bohr radius in m
  extern const double bohr_radius;
  
  //!Hartree [eV]
  extern const double Hartree;

}


#endif //_CONSTANTS_H_
