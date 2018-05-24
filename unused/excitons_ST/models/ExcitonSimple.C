// $Id: ExcitonSimple.C 4192 2015-12-10 11:11:18Z maufder $


#include "ExcitonSimple.h"
#include "SimulationInterface.h"

#include "SimulationOptions.h"
#include "Constants.h"

#include "SolveFailedException.h"

#include "Material.h"
#include "Utils.h"

#include "elem.h"
#include "ExcitonTransport.h"

#include "TiberModule.h"


ExcitonSimple::ExcitonSimple(const ModelOptions& options)
  : ExcitonProperties(options),
    _ts_r(1e45),
    _ts_nr(1e45),
    _ts_diss(1e45),
    _tt_r(1e45),
    _tt_nr(1e45),
    _tt_diss(1e45),
    _t_isc(1e45),
    _R(0.025),
    _m(1),
    _sD(1e-5),
    _tD(1e-5)
{
}


void
ExcitonSimple::do_diffusion(void)
{
  s_diffusion = _sD;
  t_diffusion = _tD;
}

void
ExcitonSimple::do_recombination()
{
  double inv_tau_s = 1.0 / _ts_nr + 1.0 / _ts_r + 1.0 / _ts_diss + 1.0 / _t_isc;
  double inv_tau_t = 1.0 / _tt_nr + 1.0 / _tt_r + 1.0 / _tt_diss;
  s_net_recomb_rate = s_density * inv_tau_s;
  t_net_recomb_rate = t_density * inv_tau_t;

  s_recombination_rate_derivative = inv_tau_s;
  t_recombination_rate_derivative = inv_tau_t;

  s_net_recomb_rate -= get_s_generation_rate();
  t_net_recomb_rate -= get_t_generation_rate() + s_density / _t_isc;

}


double ExcitonSimple::get_s_generation_rate() 
{
  std::vector<double> G(1);
  _dd_sim->get_solution(get_element(), _gen_model, G,
      std::vector<Point>(1, get_coordinates()));
  return 0.25 * G[0];
}

double ExcitonSimple::get_t_generation_rate() 
{
  std::vector<double> G(1);
  _dd_sim->get_solution(get_element(), _gen_model, G,
      std::vector<Point>(1, get_coordinates()));
  return 0.75 * G[0];
}

double
ExcitonSimple::get_s_nonradiative_recombination_rate(void)
{
  return get_s_density() / _ts_nr;
}

double
ExcitonSimple::get_t_nonradiative_recombination_rate(void)
{
  return get_t_density() / _tt_nr;
}

double
ExcitonSimple::get_s_radiative_recombination_rate(void)
{
  return get_s_density() / _ts_r ;
}

double
ExcitonSimple::get_t_radiative_recombination_rate(void)
{
  return get_t_density() / _tt_r ;
}

double
ExcitonSimple::get_s_dissociation_rate(void)
{
  return get_s_density() / _ts_diss;
}

double
ExcitonSimple::get_isc_rate(void)
{
  return get_s_density() / _t_isc;
}

double
ExcitonSimple::get_isc_rate_derivative(void)
{
  return  -1.0 / _t_isc;
}

double
ExcitonSimple::get_t_dissociation_rate(void)
{
  return get_t_density() / _tt_diss;
}


void
ExcitonSimple::prepare_element_data(void)
{

  if (!_dd_sim->is_solved())
    throw (SolveFailedException("ExcitonSimple needs solved DriftDiffusion."));

  std::vector<double> Eg(1);
  std::map<const Elem*, double>::iterator it(_bandgap_data.find(get_element()));
  if (it != _bandgap_data.end())
    set_energy(it->second);
  else
  {
    bool ok = _dd_sim->get_solution(get_element(),
        _Eg_id, Eg, std::vector<Point>(1, get_element()->centroid()));

    if (!ok)
      throw SolveFailedException("Exciton model needs everywhere drift-diffusion.");

    set_energy(Eg[0] - _R);
    _bandgap_data[get_element()] = Eg[0] - _R;
  }

  // lattice_vt was taken from TemperatureInterface
  exciton_vt = lattice_vt;

  // it's for cm
  double DOS = 3 * std::pow(2 * M_PI * Constants::me * lattice_vt * _m /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;

  set_density_of_states(DOS);
}


void
ExcitonSimple::read_database(void)
{
}


void
ExcitonSimple::do_init(void)
{
  _ts_r = get_option("tau_rad", 1e45);
  _ts_nr = get_option("tau_nonrad", 1e45);
  _ts_diss = get_option("tau_diss", 1e45);
  _sD = get_option("diffusion", 1e-5);
  _t_isc = get_option("tau_isc", 1e45);

  _tt_r = _ts_r;
  _tt_nr = _ts_nr;
  _tt_diss = _ts_diss;
  _tD = _sD;

  get_parameter("R", _R);
  get_parameter("eff_mass", _m);

  ModelOptions::submodel_iterator it_s(get_options().submodels_begin("singlet"));
  ModelOptions::submodel_iterator end_s(get_options().submodels_end("singlet"));

  if (it_s != end_s)
  {
    _ts_r = (it_s->second).get_option("tau_rad", 1e45);
    _ts_nr = (it_s->second).get_option("tau_nonrad", 1e45);
    _ts_diss = (it_s->second).get_option("tau_diss", 1e45);
    _sD = (it_s->second).get_option("diffusion", 1e-5);
  }

  ModelOptions::submodel_iterator it_t(get_options().submodels_begin("triplet"));
  ModelOptions::submodel_iterator end_t(get_options().submodels_end("triplet"));

  if (it_t != end_t)
  {
    _tt_r = (it_t->second).get_option("tau_rad", 1e45);
    _tt_nr = (it_t->second).get_option("tau_nonrad", 1e45);
    _tt_diss = (it_t->second).get_option("tau_diss", 1e45);
    _tD = (it_t->second).get_option("diffusion", 1e-5);
  }

  std::string dd = get_option("driftdiffusion_simulation", "");

  // find the drift-diffusion simulation to use
  _dd_sim = SimulationInterface::find_simulation(dd);

  if (_dd_sim == NULL)
  {
    std::string msg("ExcitonSimple: Simulation " +
        std::string(dd) + " not found");
    throw InitFailedException(msg);
  }

  _gen_model = _dd_sim->get_solution_id("eExcitonGeneration");

  _Eg_id = _dd_sim->get_solution_id("Eg");

  if ((_gen_model == INVALID_ID) || (_Eg_id == INVALID_ID))
  {
    std::string msg("ExcitonSimple: Simulation " +
        std::string(dd) + " does not provide all necessary variables");
    throw InitFailedException(msg);
  }

}

