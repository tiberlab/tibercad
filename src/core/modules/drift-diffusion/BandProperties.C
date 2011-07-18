// $Id$

#include "BandProperties.h"
#include "DriftDiffusionProperties.h"
#include "ModelOptions.h"

const double
BandProperties::_dos_factor = pow(2.0 * M_PI *
    Constants::me / (Constants::h * Constants::h) *
    Constants::e, 1.5) / 1e6;

BandProperties::BandProperties() : DriftDiffusionModelInterface(ModelOptions()) {}


BandProperties::BandProperties(const ModelOptions& options) :
    DriftDiffusionModelInterface(options)
{
  has_option("particle");
}

double
BandProperties::get_lattice_temperature(void) const
{
  return get_driftdiffusionproperties().get_lattice_temperature();
}


std::pair<double, double>
BandProperties::get_density_and_derivative(void) const
{
  return std::pair<double, double>(0,0);
}
