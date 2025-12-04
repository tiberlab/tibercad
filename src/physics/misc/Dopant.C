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
 * \file Dopant.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/misc/Dopant.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/profiles/ExternalProfile.h"

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

  dop = new Dopant(density, ionisation_energy, g_factor, type);
  dop->_options = options;
  dop->_incomplete_ionisation = options.get_option("incomplete_ionization", true);

  if (options.has_submodel("profile"))
  {
    dop->_profile = ExternalProfile::create(
        options.submodels_begin("profile")->second);
  }

  return dop;
}


void
Dopant::calculate_doping_density(const libMesh::Elem* elem, const libMesh::Point& p)
{
  if (_profile != nullptr)
    _density = _profile->get_data(elem, p);
}




double
Dopant::get_ionized_dopant_density(double arg, double kT)
{
  if (_incomplete_ionisation)
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
  if (_incomplete_ionisation)
  {
    double tmp = _g_factor * std::exp((arg + _ionisation_energy) / kT);
    double denom = 1.0 + tmp;
    denom *= denom;

    return -tmp * _density / (denom * kT);
  }
  else
    return 0.0;
}



