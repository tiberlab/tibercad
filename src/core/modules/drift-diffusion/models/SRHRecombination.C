// $Id$

#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"



TIBER_MODULE(SRHRecombination, recombination, srh)



void
SRHRecombination::read_database(void)
{
  Database& db = get_database();

  if (get_option("trap", false))
    db.set_section("recombination/trap");
  else
    db.set_section("recombination/SRH");

  _E_t = db.get("Etrap", _E_t);

  std::vector<double> data(2, 0);

  db.get("Talpha", data);
  _Talpha_e = data[0];
  _Talpha_h = data[1];

  data = std::vector<double>(2, 0);
  db.get("Tcoeff", data);
  _Tcoeff_e = data[0];
  _Tcoeff_h = data[1];

  data[0] = data[1] = 1e-12;
  db.get("taumin", data);
  double taumin_e = data[0];
  double taumin_h = data[1];
  data[0] = data[1] = 1e-9;
  db.get("taumax", data);
  double taumax_e = data[0];
  double taumax_h = data[1];
  data[0] = data[1] = 1e16;
  db.get("Nref", data);
  double Nref_e = data[0];
  double Nref_h = data[1];
  data[0] = data[1] = 2;
  db.get("gamma", data);
  double g_e = data[0];
  double g_h = data[1];

  double N = get_material()->get_total_doping_density();

  // electrons
  double denom = 1.0 + std::pow(N / Nref_e, g_e);
  _tau_n = taumin_e + (taumax_e - taumin_e) / denom;

  // holes
  denom = 1.0 + std::pow(N / Nref_h, g_h);
  _tau_p = taumin_h + (taumax_h - taumin_h) / denom;
}




void
SRHRecombination::read_interface_database(void)
{
  Database& db = get_database();

  db.set_section("recombination/surface_rec");

  _E_t = db.get("Etrap", _E_t);

  std::vector<double> data(2, 1e7);
  db.get("rec_velocity", data);

  _tau_n = 1.0 / data[0];
  _tau_p = 1.0 / data[1];

}


void
SRHRecombination::do_init(void)
{
  if (get_option("trap", false))
  {
    get_parameter("sigma_n", _sigma_n);
    get_parameter("sigma_p", _sigma_p);
    _trap = true;
  }

  get_parameter("tau_n", _tau_n);
  get_parameter("tau_p", _tau_p);

  get_parameter("Et", _E_t);
  get_parameter("Nt", _density);

  std::string tmp("m");
  tmp[0] = _energy_reference;
  tmp =  get_option("reference", tmp);
  _energy_reference = tmp[0];

}



void
SRHRecombination::do_init_interface(const Material* comp_A,
    const Material* comp_B)
{
  get_parameter("rec_velocity_n", _tau_n, true, new Invert());
  get_parameter("rec_velocity_p", _tau_p, true, new Invert());

  do_init();

}


double
SRHRecombination::get_trap_level(void)
{
  double ref;
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Ec = dd.get_conduction_band_edge();
  double Ev = dd.get_valence_band_edge();
  switch (_energy_reference)
  {
    case 'v':
      ref = Ev + _E_t;
      break;

    case 'm':
      ref = 0.5 * (Ev + Ec) + _E_t;
      break;

    default:
      ref = Ec - _E_t;
      break;
  }

  return ref;
}



void
SRHRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  double T = dd.get_lattice_temperature();
  // TODO should take carrier temperatures

  double Et = get_trap_level() - dd.get_equilibrium_fermi_level();
  double f = std::exp(Et / T);

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    tau_n = 1.0 / (_density * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (_density * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  double denom = tau_p * (n + gn * n0 * f) + tau_n * (p + gp * p0 / f);
  double tmp = n * p / denom;
  recomb_e = recomb_h = tmp - gn * gp * n0 * p0 / denom;
}



void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  double T = dd.get_lattice_temperature();

  double Et = get_trap_level() - dd.get_equilibrium_fermi_level();
  double f = std::exp(Et / T);

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    tau_n = 1.0 / (_density * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (_density * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  long double denom = tau_p * (n + gn * n0 * f) + tau_n * (p + gp * p0 / f);
  long double tmp = n * p / denom;
  long double SRH = tmp - gn * gp * n0 * p0 / denom;

  long double a = p / denom;
  a = a - tau_p * SRH / denom;
  long double b = n / denom;
  b = b - tau_n * SRH / denom;

  if (a < 0) a = p / denom;
  if (b < 0) b = n / denom;

  recomb_e[0] = recomb_h[0] = a;
  recomb_e[1] = recomb_h[1] = b;
}


void
SRHRecombination::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const SRHRecombination* scA =
    dynamic_cast<const SRHRecombination*>(comp_A);
  const SRHRecombination* scB =
    dynamic_cast<const SRHRecombination*>(comp_B);

  _tau_n = alloy(scA->_tau_n, scB->_tau_n, xa);
  _tau_p = alloy(scA->_tau_p, scB->_tau_p, xa);
}

