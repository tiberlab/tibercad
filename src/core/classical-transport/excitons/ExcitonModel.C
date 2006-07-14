// $Id$


#include "ExcitonModel.h"
#include "RecombinationModelInterface.h"

#include "SemiconductorModel.h"
#include "SimulationOptions.h"
#include "Constants.h"
#include <iostream>

ExcitonModel::ExcitonModel(void)
  : ExcitonProperties()
{
}


void
ExcitonModel::calculate_all(double fermi_x, const Point& coord)
{
  double kT = exciton_vt;

  // 1.) exciton density
  double arg_x = - (band_gap - _R - fermi_x) / kT;

  calculate_density_and_derivative (arg_x, density, density_derivative);

  double Nx = _DOS;
  density *= Nx;
  density_derivative *= Nx / kT;


  // 2.) mobility
  //   mu_x = e * D_x /(kB*T) = D_x/kT
  mobility = _mu;

  // 3.) Recombination
  recombination_rate = density / _t;
  recombination_rate_derivative = density_derivative / _t;

  // 4.) Generation
  RecombinationModelInterface* rec =
    get_driftdiffusion_properties()->get_recombination_model(_gen_mod_id);
  assert(rec != NULL);
  // the Drift-Diffusion recombination models return a net recombination rate
  // for electrons and holes.
  double dummy;
  rec->get_net_recombination_rates(generation_rate, dummy);

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
  // find out which ID has the exciton generation model
  // for now the model name is hard wired
  _gen_mod_id = RecombinationModelInterface::get_id(_exciton_generation_model);

  // TODO perhaps we should check if the exciton generation model exists?
  assert(_gen_mod_id != 0);
}
