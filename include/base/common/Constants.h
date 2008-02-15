// $Id$

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

//! Physical constants
/*!
 *
 * Physical constants are given in SI units
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

  //!Lorentz Numner [(W * Ohm) / ( K * K) ]
  extern const double Lorenz_Number;


  //!Speed of light [m/s]
  extern const double c;


  //!Atomic unit of time [second] \f$ \hbar / Ha \f$
  extern const double atomic_time;


  //!Fine structure constant. In Gauss units: \f$ e^2/(\hbar c)  \f$ In SI units: \f$  e^2/(\hbar c 4 \pi \varepsilon_0)  \f$
  extern const double fine_structure_constant;

 
 
  //!Unit of field c/1e6  (300B/cm)
  extern const double field_gauss_unit;

}


#endif //_CONSTANTS_H_
