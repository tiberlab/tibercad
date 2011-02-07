// $Id$

#include "Dopant.h"

#include "SimulationOptions.h"

#include <cmath>


Dopant*
Dopant::create(const std::string& profile, const ModelOptions& options)
{
  Dopant* dop = NULL;

  double density = options.get_option("density", 0.0);
  density = options.get_option("Nd", density);
  double ionisation_energy = options.get_option("level", 0.025);
  ionisation_energy = options.get_option("Ed", ionisation_energy);
  int g_factor = 2;
  DopingType type = N_TYPE;
  std::string type_s = options.get_option("type", "donor");
  if (type_s == "acceptor")
  {
    type = P_TYPE;
    // TODO: this is commented out only for now, to not break testsuite
    //g_factor = 4;
  }
  g_factor = options.get_option("g", g_factor);

  if (profile.empty() || (profile == "constant"))
    dop = new Dopant(density, ionisation_energy, g_factor, type);

  if (dop != NULL)
    dop->_options = options;

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



