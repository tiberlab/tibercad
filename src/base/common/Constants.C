// $Id$

#include "tibercad/physics/Constants.h"

using namespace std;

const double
Constants::k_Boltzmann = 8.617343e-5;

const double&
Constants::k_B = k_Boltzmann;

const double&
Constants::kb = k_Boltzmann;

const double
Constants::elementary_charge = 1.602176487e-19;

const double&
Constants::e = elementary_charge;

const double
Constants::electron_volt = elementary_charge;

const double&
Constants::eV = electron_volt;

const double
Constants::epsilon = 8.854187871e-12;

const double&
Constants::eps = epsilon;

const double&
Constants::e0 = epsilon;

const double
Constants::electron_mass = 9.10938215e-31;

const double&
Constants::me = electron_mass;

const double
Constants::plancks_constant = 6.62606896e-34;

const double&
Constants::h = plancks_constant;

const double
Constants::hbar = h / (2 * M_PI); // 1.0545716e-34

const double
Constants::bohr_radius = 4 * M_PI * epsilon * (hbar * hbar) /
  (electron_mass * elementary_charge * elementary_charge); // 5.2917721e-11

const double
Constants::Hartree = (hbar * hbar) /
  (electron_mass * bohr_radius * bohr_radius) / elementary_charge; // 27.2114

const double
Constants::Lorenz_Number = (4.0 * k_Boltzmann *  k_Boltzmann ) / (M_PI);

const double
Constants::c = 299792458.0; 


const double
Constants::atomic_time = hbar/(Hartree*elementary_charge);

const double
Constants::fine_structure_constant = (elementary_charge * elementary_charge) 
  /( plancks_constant/(2.0*M_PI) * c * 4.0 * M_PI * epsilon ); 



const double
Constants::field_gauss_unit = 1/(1e6 /(c * 1e2));


const double
Constants::avogadro = 6.02214179e23;

const double
Constants::EE = 1.439966289;

