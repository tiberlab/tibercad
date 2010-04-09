// $Id$


#include "ExcitonModel.h"
#include "SimulationInterface.h"

#include "SimulationOptions.h"
#include "Constants.h"

#include "SolveFailedException.h"

#include "Material.h"
#include "Utils.h"

#include "elem.h"


TIBER_MODULE(ExcitonModel, exbulk, simple)



ExcitonModel::ExcitonModel(const ModelOptions& options)
  : ExcitonProperties(options),
    _t_r(1e45),
    _t_nr(1e45),
    _t_diss(1e45),
    _R(0.025),
    _m(1),
    _mu(1500)
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


  double G;
  _dd_sim->get_solution(get_element(), get_coordinates(), _gen_model, G);
  net_recomb_rate -= G;
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

  double Eg;
  std::map<const Elem*, double>::iterator it(_bandgap_data.find(get_element()));
  if (it != _bandgap_data.end())
    set_energy(it->second);
  else
  {
    bool ok = ((SimulationInterface*)_dd_sim)->get_solution(get_element(),
        get_element()->centroid(), _Eg_id, Eg);

    if (!ok)
      throw SolveFailedException("Exciton model needs everywhere drift-diffusion.");

    set_energy(Eg - _R);
    _bandgap_data[get_element()] = Eg - _R;
  }

  // lattice_vt was taken from TemperatureInterface
  exciton_vt = lattice_vt;

  // it's for cm
  double DOS = 3 * std::pow(2 * M_PI * Constants::me * lattice_vt * _m /
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
  get_parameter("tau_rad", _t_r);
  get_parameter("tau_nonrad", _t_nr);
  get_parameter("tau_diss", _t_diss);

  get_parameter("R", _R);
  get_parameter("eff_mass", _m);
  get_parameter("mobility", _mu);


  std::string dd = get_option("DD_simulation", "driftdiffusion");

  // find the drift-diffusion simulation to use
  _dd_sim = SimulationInterface::find_simulation(dd);

  if (_dd_sim == NULL)
  {
    std::string msg("ExcitonModel: Simulation " +
        std::string(dd) + " not found");
    throw InitFailedException(msg);
  }

  std::string varname("recombination.");
  varname += get_option("generation_model", "");
  _gen_model = _dd_sim->get_solution_id(varname);

  _Eg_id = _dd_sim->get_solution_id("Eg");

  if ((_gen_model == INVALID_ID) || (_Eg_id == INVALID_ID))
  {
    std::string msg("ExcitonModel: Simulation " +
        std::string(dd) + " does not provide all necessary variables");
    throw InitFailedException(msg);
  }

}

