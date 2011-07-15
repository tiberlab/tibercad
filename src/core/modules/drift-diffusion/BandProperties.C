// $Id$

#include "BandProperties.h"
#include "ModelOptions.h"

const double
BandProperties::_dos_factor = pow(2.0 * M_PI *
    Constants::me / (Constants::h * Constants::h) *
    Constants::e, 1.5) / 1e6;

BandProperties::BandProperties() : DriftDiffusionModelInterface(ModelOptions()) {}


BandProperties::BandProperties(const ModelOptions& options) :
    DriftDiffusionModelInterface(options)
{

}

