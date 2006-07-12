#include "Dopant.h"

#include <cmath>

double
Dopant::get_ionized_dopant_density(double arg, double kT)
{
  double denom = 1.0 + _g_factor * std::exp((arg + _ionisation_energy) / kT);

  return _density / denom;
}


double
Dopant::get_ionized_dopant_density_derivative(double arg, double kT)
{
  double tmp = _g_factor * std::exp((arg + _ionisation_energy) / kT);
  double denom = 1.0 + tmp;
  denom *= denom;

  return -tmp * _density / (denom * kT);
}



