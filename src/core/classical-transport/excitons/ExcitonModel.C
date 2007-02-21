// $Id$


#include "ExcitonModel.h"
#include "RecombinationModelInterface.h"

#include "SemiconductorModel.h"
#include "SimulationOptions.h"
#include "Constants.h"

#include "SolveFailedException.h"

#include "DriftDiffusion.h"
#include "Material.h"
#include "Utils.h"


ExcitonModel::ExcitonModel(void)
{
}


void
ExcitonModel::do_mobility(void)
{
  mobility = _mu;
}

void
ExcitonModel::do_recombination(void)
{
  double inv_tau = 1.0 / _t_nr + 1.0 / _t_r + 1.0 / _t_diss;
  net_recomb_rate = density * inv_tau;
  recombination_rate_derivative = density_derivative * inv_tau;

  ID dd_id = _dd_sim->get_id();
  DriftDiffusionProperties* ddprop =
    static_cast<DriftDiffusionProperties*>(get_material()->get_model(dd_id));

  
  double kT = ddprop->get_lattice_temperature();
  ddprop->set_carrier_temperatures(kT, kT);

  DriftDiffusion::Solution sol;
  _dd_sim->get_solution(get_element(), get_coordinates(), sol);
  ddprop->set_potentials(sol.potential, sol.fermi_e, sol.fermi_h);
  ddprop->calculate_densities();
  RecombinationModelInterface* recmod =
    ddprop->get_recombination_model(_gen_model);
  if (recmod != NULL)
  {
    double G, dummy;
    recmod->get_net_recombination_rates(G, dummy);
    net_recomb_rate -= G;
  }
}


double
ExcitonModel::get_nonradiative_recombination_rate(void)
{
  return density / _t_nr;
}


double
ExcitonModel::get_radiative_recombination_rate(void)
{
  return density / _t_r;
}


double
ExcitonModel::get_dissociation_rate(void)
{
  return density / _t_diss;
}


void
ExcitonModel::prepare_element_data(void)
{

  if (!_dd_sim->is_solved())
    throw (SolveFailedException("ExcitonModel needs solved DriftDiffusion."));

  // reinit DriftDiffusionProperties
  ID dd_id = _dd_sim->get_id();
  DriftDiffusionProperties* ddprop =
    static_cast<DriftDiffusionProperties*>(get_material()->get_model(dd_id));
  ddprop->reinit(get_element());

  // this does not depend on the coordinate
  set_energy(ddprop->get_band_gap() - _R);

  double kT = ddprop->get_lattice_temperature();
  exciton_vt = lattice_vt = kT;

  // it's for cm
  double DOS = 3 * std::pow(2 * M_PI * Constants::me * kT * _m /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;

  set_density_of_states(DOS);
}


void
ExcitonModel::read_database(void)
{
}


void
ExcitonModel::do_init(void)
{

  _t_r = get_options().get_option("tau_rad", 1.0e45);
  _t_nr = get_options().get_option("tau_nonrad", 1.0e45);
  _t_diss = get_options().get_option("tau_diss", 1.0e45);

  _R = get_options().get_option("R", 0.025);
  _m = get_options().get_option("eff_mass", 1.0);
  _mu = get_options().get_option("mobility", 1500.0);

  // look for the generation model
  // (default is no generation at all)
  _gen_model = get_id_from_name<RecombinationModelInterface>(
    get_options().get_option("generation_model", ""));

  std::string dd = get_options().get_option("DD_simulation",
      Utils::extract_typename(typeid(_dd_sim)));

  // find the drift-diffusion simulation to use
  _dd_sim = dynamic_cast<DriftDiffusion*>(
      SimulationInterface::find_simulation(dd));

  if (_dd_sim == NULL)
  {
    std::string msg("ExcitonModel: Simulation " +
        std::string(dd) + " not found");
    throw InitFailedException(msg);
  }

}


void
ExcitonModel::copy_from(const PhysicalModelInterface* rhs)
{
  ExcitonProperties::copy_from(rhs);

  const ExcitonModel* mod =
    dynamic_cast<const ExcitonModel*>(rhs);

  _t_r = mod->_t_r;
  _t_nr = mod->_t_nr;
  _t_diss = mod->_t_diss;
  _R = mod->_R;
  _m = mod->_m;
  _mu = mod->_mu;
  _gen_model = mod->_gen_model;
  _dd_sim = mod->_dd_sim;
}


