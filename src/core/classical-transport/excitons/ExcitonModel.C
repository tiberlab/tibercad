// $Id$


#include "ExcitonModel.h"

#include "SemiconductorModel.h"
#include "SimulationOptions.h"
#include "Constants.h"

ExcitonModel::ExcitonModel(void)
  : ExcitonProperties()
{
}

void
ExcitonModel::calculate_all(double fermi_x, const Point& coord)
{
  double kT = exciton_vt;

}

void
ExcitonModel::prepare_element_data(void)
{
  double kT = SimulationOptions::T * Constants::k_B;
  exciton_vt = kT;
  // it's for cm
  _DOS = 3 * std::pow(2 * M_PI * Constants::me * kT * _m /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;

  band_gap = get_driftdiffusion_properties()->get_band_gap();
}

void
ExcitonModel::read_database(const Dummy& d)
{
}
