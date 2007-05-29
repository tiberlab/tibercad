// $Id$

#include "Dopant.h"

#include "SimulationOptions.h"

#include <cmath>

double
Dopant::get_ionized_dopant_density(double arg, double kT)
{
  if (SimulationOptions::incomplete_ionization)
  {
    double denom = 1.0 + _g_factor * std::exp((arg + _ionisation_energy) / kT);
    return _density / denom;
  }
  else
    return _density;
}


double
Dopant::get_ionized_dopant_density_derivative(double arg, double kT)
{
  if (SimulationOptions::incomplete_ionization)
  {
    double tmp = _g_factor * std::exp((arg + _ionisation_energy) / kT);
    double denom = 1.0 + tmp;
    denom *= denom;

    return -tmp * _density / (denom * kT);
  }
  else
    return 0.0;
}



