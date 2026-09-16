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
 * \file Constants.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_CONSTANTS_H
#define TC_CONSTANTS_H

#include "tibercad/base/tiber_dll.h"

#ifndef M_PI
 #define M_PI		3.14159265358979323846
#endif
#ifndef M_2_SQRTPI	
 #define M_2_SQRTPI	1.12837916709551257390
#endif

//! Physical constants
/*!
 *
 * Physical constants are given in SI units
 *
 */
namespace Constants 
{
  
  //! Boltzmann constant in eV / K
  extern const double k_Boltzmann;
  extern const double& k_B;
  extern const double& kb;

  //! Elementary charge in C
  extern const double elementary_charge;
  extern const double& e;

  //! 1 eV
  extern const double electron_volt;
  extern const double& eV;

  //! dielectric constant of vacuum
  extern const double epsilon;
  extern const double& eps;
  extern const double& e0;

  //! rest mass of an electron in kg
  extern const double electron_mass;
  extern const double& me;

  //! Planck's constant in eV * s
  extern const double plancks_constant;
  extern const double& h;
  extern const double hbar;


  
  //! Bohr radius in m
  extern const double bohr_radius;
  
  //! Hartree [eV]
  extern const double Hartree;

  //! Lorentz Numner [(W * Ohm) / ( K * K) ]
  extern const double Lorenz_Number;


  //! Speed of light [m/s]
  extern const double c;


  //! Atomic unit of time [second] \f$ \hbar / Ha \f$ (is about 2.4e-17)
  extern const double atomic_time;


  //! Fine structure constant. In Gauss units: \f$ e^2/(\hbar c)  \f$ In SI units: \f$  e^2/(\hbar c 4 \pi \varepsilon_0)  \f$
  extern const double fine_structure_constant;

 
 
  //!Unit of field c/1e6  (300B/cm)
  extern const double field_gauss_unit;



  //! Avogadro's number
  extern const double avogadro;

  //! e^2/(4 pi e0) in eV * nm  
  extern const double EE;
}


#endif //_CONSTANTS_H_
