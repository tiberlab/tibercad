// $Id$

#include "Constants.h"
#include <cmath>

const double
Constants::k_Boltzmann = 8.617343e-5;

const double&
Constants::k_B = k_Boltzmann;

const double&
Constants::kb = k_Boltzmann;

const double
Constants::elementary_charge = 1.60219e-19;

const double&
Constants::e = elementary_charge;

const double
Constants::electron_volt = elementary_charge;

const double&
Constants::eV = electron_volt;

const double
Constants::epsilon = 8.85374e-12;

const double&
Constants::eps = epsilon;

const double&
Constants::e0 = epsilon;

const double
Constants::electron_mass = 9.10956e-31;

const double&
Constants::me = electron_mass;

const double
Constants::plancks_constant = 6.62620e-34;

const double&
Constants::h = plancks_constant;

const double
Constants::hbar = h / (2 * M_PI);

const double
Constants::bohr_radius = 4 * M_PI * epsilon * (hbar * hbar) /(electron_mass * elementary_charge * elementary_charge);

const double 
Constants::Hartree = (hbar * hbar)/(electron_mass * bohr_radius * bohr_radius) / elementary_charge;
