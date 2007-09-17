// $Id$

#include "Dopant.h"

#include "SimulationOptions.h"

#include <cmath>


Dopant*
Dopant::create(const std::string& profile, const ModelOptions& options)
{
  Dopant* dop = NULL;

  double density = options.get_option("density", 0.0);
  double ionisation_energy = options.get_option("level", 0.025);
  int g_factor = options.get_option("g", 2);
  std::string type_s = options.get_option("type", "donor");
  DopingType type = N_TYPE;
  if (type_s == "acceptor")
    type = P_TYPE;

  if (profile == "constant")
    dop = new Dopant(density, ionisation_energy, g_factor, type);

  if (dop != NULL)
  {
    dop->_options = options;
    (dop->_options).delete_option("density");
    (dop->_options).delete_option("level");
    (dop->_options).delete_option("g");
    (dop->_options).delete_option("type");
  }

  return dop;
}


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



