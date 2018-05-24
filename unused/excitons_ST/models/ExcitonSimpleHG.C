#include "ExcitonSimpleHG.h"
#include "SimulationInterface.h"

#include "SimulationOptions.h"
#include "Constants.h"
#include "Database.h"
#include "SolveFailedException.h"

#include "Material.h"
#include "Utils.h"

#include "elem.h"
#include "ExcitonTransport.h"

#include "TiberModule.h"

bool ExcitonSimpleHG::_coupled;

ExcitonSimpleHG::ExcitonSimpleHG(const ModelOptions& options)
  : ExcitonProperties(options),
    _ts_r(1e-9),
    _ts_nr(1e-9),
    _ts_diss(1e-9),
    _tt_r(1e-9),
    _tt_nr(1e-9),
    _tt_diss(1e-9),
    _t_isc(1e45),
    _R(0.025),
    _m(1.0),
    _sD(1e-5),
    _tD(1e-5),
    _Rf(3e-7),
    _Rd(1e-7),
    _R_hg(1.5e-7)
{
}


void
ExcitonSimpleHG::do_diffusion(void)
{
  s_diffusion = _sD;
  t_diffusion = _tD;
}

void
ExcitonSimpleHG::do_recombination()
{
  double inv_tau_s = 1.0 / _ts_nr + 1.0 / _ts_r + 1.0 / _ts_diss + 1.0 / _t_isc;
  double inv_tau_t = 1.0 / _tt_nr + 1.0 / _tt_r + 1.0 / _tt_diss;

  if (_hg_sim != NULL)
  {
    if (_hg_sim->is_solved()) 
    {
      inv_tau_s += _Kf + _Kds;
      inv_tau_t += _Kdt;
    }
  }

  s_net_recomb_rate = s_density * inv_tau_s;
  t_net_recomb_rate = t_density * inv_tau_t;

  s_recombination_rate_derivative = inv_tau_s;
  t_recombination_rate_derivative = inv_tau_t;

  s_net_recomb_rate -= get_s_generation_rate() + get_s_hg_generation_rate();
  t_net_recomb_rate -= get_t_generation_rate() + get_t_hg_generation_rate() + s_density / _t_isc;

}


double ExcitonSimpleHG::get_s_generation_rate() 
{
  std::vector<double> G(1);
  _dd_sim->get_solution(get_element(), _gen_model, G,
      std::vector<Point>(1, get_coordinates()));
  return 0.25 * G[0];
}

double ExcitonSimpleHG::get_s_hg_generation_rate() 
{
  if (_hg_sim != NULL)
  {
    if (_hg_sim->is_solved() && _coupled) 
    {
      std::vector<double> G(1);
      _coupled = false;
      _hg_sim->get_solution(get_element(), _exs_model, G, std::vector<Point>(1, get_coordinates()));
      _coupled = true;
      return G[0];
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}

double ExcitonSimpleHG::get_t_generation_rate() 
{
  std::vector<double> G(1);
  _dd_sim->get_solution(get_element(), _gen_model, G, std::vector<Point>(1, get_coordinates()));
  return 0.75 * G[0];
}

double ExcitonSimpleHG::get_t_hg_generation_rate() 
{
  if (_hg_sim != NULL)
  {
    if (_hg_sim->is_solved() && _coupled) 
    {
      std::vector<double> G(1);
      _coupled = false;
      _hg_sim->get_solution(get_element(), _ext_model, G, std::vector<Point>(1, get_coordinates()));
      _coupled = true;
      return G[0];
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}

double
ExcitonSimpleHG::get_s_nonradiative_recombination_rate(void)
{
  return get_s_density() / _ts_nr;
}

double
ExcitonSimpleHG::get_t_nonradiative_recombination_rate(void)
{
  return get_t_density() / _tt_nr;
}

double
ExcitonSimpleHG::get_s_radiative_recombination_rate(void)
{
  return get_s_density() / _ts_r ;
}

double
ExcitonSimpleHG::get_t_radiative_recombination_rate(void)
{
  return get_t_density() / _tt_r ;
}

double
ExcitonSimpleHG::get_s_dissociation_rate(void)
{
  return get_s_density() / _ts_diss;
}

double
ExcitonSimpleHG::get_s_hg_recombination_rate(void)
{
  if (_hg_sim != NULL)
  {
    if (_hg_sim->is_solved()) 
    {
      return get_s_density() * (_Kf + _Kds);
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}

double
ExcitonSimpleHG::get_t_hg_recombination_rate(void)
{
  if (_hg_sim != NULL)
  {
    if (_hg_sim->is_solved()) 
    {
      return get_t_density() * _Kdt;
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;  
  }
}

double
ExcitonSimpleHG::get_isc_rate(void)
{
  return get_s_density() / _t_isc;
}

double
ExcitonSimpleHG::get_isc_rate_derivative(void)
{
  return  -1.0 / _t_isc;
}

double
ExcitonSimpleHG::get_t_dissociation_rate(void)
{
  return get_t_density() / _tt_diss;
}


void
ExcitonSimpleHG::prepare_element_data(void)
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
ExcitonSimpleHG::read_database(void)
{
  const Database& db = get_database();

  db.set_section("permittivity");
  _er = db.get("permittivity", 1.0);

}


void
ExcitonSimpleHG::do_init(void)
{
  _ts_r = get_option("tau_rad", 1e-9);
  _ts_nr = get_option("tau_nonrad", 1e-9);
  _ts_diss = get_option("tau_diss", 1e-9);
  _sD = get_option("diffusion", 1e-5);
  _t_isc = get_option("tau_isc", 1e45);

  _Rf   = get_option("foster_radius", 3e-7);
  _Rd   = get_option("dexter_radius", 1e-7);
  _R_hg = get_option("average_distance", 1.5e-7);

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
    _ts_r = (it_s->second).get_option("tau_rad", 1e-9);
    _ts_nr = (it_s->second).get_option("tau_nonrad", 1e-9);
    _ts_diss = (it_s->second).get_option("tau_diss", 1e-9);
    _sD = (it_s->second).get_option("diffusion", 1e-5);
  }

  ModelOptions::submodel_iterator it_t(get_options().submodels_begin("triplet"));
  ModelOptions::submodel_iterator end_t(get_options().submodels_end("triplet"));

  if (it_t != end_t)
  {
    _tt_r = (it_t->second).get_option("tau_rad", 1e-9);
    _tt_nr = (it_t->second).get_option("tau_nonrad", 1e-9);
    _tt_diss = (it_t->second).get_option("tau_diss", 1e-9);
    _tD = (it_t->second).get_option("diffusion", 1e-5);
  }

  std::string dd = get_option("driftdiffusion_simulation", "");
  std::string hg = get_option("hostguest_simulation", "");

  // find the drift-diffusion simulation to use
  _dd_sim = SimulationInterface::find_simulation(dd);
  _hg_sim = SimulationInterface::find_simulation(hg);

  if (_dd_sim == NULL)
  {
    std::string msg("ExcitonSimpleHG: Simulation " +
        std::string(dd) + " not found");
    throw InitFailedException(msg);
  }

  _gen_model = _dd_sim->get_solution_id("eExcitonGeneration");

  _Eg_id = _dd_sim->get_solution_id("Eg");

  if ((_gen_model == INVALID_ID) || (_Eg_id == INVALID_ID))
  {
    std::string msg("ExcitonSimpleHG: Simulation " +
        std::string(dd) + " does not provide all necessary variables");
    throw InitFailedException(msg);
  }

  _coupled = false;
  
  if (_hg_sim == NULL)
  {
    std::string msg("ExcitonSimpleHG: host-guest simulation " +
        std::string(hg) + " not found");
    throw InitFailedException(msg);
  }
  

  if (_hg_sim != NULL)
  {
    _exs_model = _hg_sim->get_solution_id("s_hg_recombination");
    _ext_model = _hg_sim->get_solution_id("t_hg_recombination");

    _coupled = true;
  }

  double RBeff = Constants::bohr_radius * 100.0 * _er / _m;  //effective Bohr radius in cm
  double inv_tau_s = 1.0 / _ts_nr + 1.0 / _ts_r + 1.0 / _ts_diss + 1.0 / _t_isc;
  double inv_tau_t = 1.0 / _tt_nr + 1.0 / _tt_r + 1.0 / _tt_diss;

  _Kf = pow(_Rf / _R_hg, 6) / _ts_r;  //Foster rate
  _Kds = _Kdt = 1e3 * exp(2*(_Rd - _R_hg) / RBeff);  //Dexter rate
  _Kds *= inv_tau_s;
  _Kdt *= inv_tau_t;

  //std::cout<<"_Kf = "<<_Kf<<" _Kds = "<<_Kds<<" _Kdt = "<< _Kdt<<" Kisc = "<<1.0 / _t_isc<<std::endl;
}

