// $Id$

#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"
#include "DensityOfStates.h"
#include "ExternalProfile.h"
#include "TiberMath.h"

#include "Material.h"
#include "Database.h"


#include "TiberModule.h"


using namespace std;

SRHRecombination::TrapAssisted::TrapAssisted(void) :
  m_trap(0.25)
{
}


double
SRHRecombination::TrapAssisted::get_gamma(double F, double T, double Et)
{
  if (Et == 0.0) return 0.0;

  double Et_kT = Et / T;
  double Kn_ref = 2.0/3.0 * Et_kT;
  double Kn = 4.0/3.0 * sqrt(2 * m_trap * Constants::electron_mass *
      Constants::e * Et * Et * Et) /
      (3 * Constants::hbar * F);

  double gamma = 0;

  if (Kn > Kn_ref)
  {
    double fg = sqrt(24 * Constants::electron_mass * m_trap * Constants::e * T * T * T) / Constants::hbar;
    gamma = 2 * sqrt(3*M_PI) * F / fg * exp(F*F/(fg*fg));
  }
  else
  {
    const double a1 = 0.3480242;
    const double a2 = -0.0958798;
    const double a3 = 0.7478556;
    //const double p1 = 0.47047; // ?
    const double p1 = 0.61685; // ?

    double ttn = 1.0 / (1 + p1 * (0.5 * Et_kT - Kn*0.75) / sqrt(Kn * 0.375));
    gamma = Et_kT * sqrt(2 * M_PI / (3 * Kn)) *
        (a1 * ttn + a2 * ttn * ttn + a3 * ttn * ttn * ttn) * exp(Et_kT - Kn);
  }

  return gamma;
}



SRHRecombination::SRHRecombination(const ModelOptions& options) :
  RecombinationModelInterface(options),
  _trap(false),
  _tau_n(1e-9),
  _tau_p(1e-9),
  _sigma_n(1e-15),
  _sigma_p(1e-15),
  _gen_TC(0.0),
  _gen_VT(0.0),
  _E_t(0.0),
  _density(1e16),
  _energy_reference('m'),
  _Talpha_e(0.0),
  _Talpha_h(0.0),
  _Tcoeff_e(0.0),
  _Tcoeff_h(0.0),
  _tat(nullptr),
  _dos(nullptr),
  _profile(nullptr)
{
}



SRHRecombination::~SRHRecombination(void)
{
  delete _tat;
  delete _profile;
  destroy(_dos);
}



void
SRHRecombination::prepare_submodels(void)
{
  if (get_options().has_submodel("density_of_states"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("density_of_states"));
    _dos = DensityOfStates::create(it->second);
    add_submodel("dos", _dos);
  }
}



void
SRHRecombination::read_database(void)
{
  const Database& db = get_database();

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
  const Database& db = get_database();

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
    if (get_options().has_submodel("profile"))
    {
      _profile = ExternalProfile::create(
          get_options().submodels_begin("profile")->second);
      _density = 1.0;
    }

    get_parameter("sigma_n", _sigma_n);
    get_parameter("sigma_p", _sigma_p);
    get_parameter("Nt", _density);

    get_parameter("trap_to_cb_rate", _gen_TC);
    get_parameter("vb_to_trap_rate", _gen_VT);

    _trap = true;
  }


  get_parameter("tau_n", _tau_n);
  get_parameter("tau_p", _tau_p);

  get_parameter("Et", _E_t);

  std::string tmp("m");
  tmp[0] = _energy_reference;
  tmp =  get_option("reference", tmp);
  _energy_reference = tmp[0];

  ModelOptions::submodel_iterator it(get_options().submodels_begin("trap_assisted_tunneling"));
  if (get_option("trap_assisted_tunneling", false) ||
      it != get_options().submodels_end("trap_assisted_tunneling"))
  {
    _tat = new TrapAssisted();

    if (it != get_options().submodels_end("trap_assisted_tunneling"))
      _tat->m_trap = it->second.get_option("tunneling_mass", _tat->m_trap);
  }
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

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double T = dd.get_lattice_temperature();

  double Et = get_trap_level();

  double dens = _density;
  if (_profile != nullptr)
    dens *= _profile->get_data(dd.get_element(), dd.get_coordinates());

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    dens = max(dens, 1e-32);
    tau_n = 1.0 / (dens * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (dens * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  if (_tat != NULL)
  {
    // definition of integration limit
    double Ec = dd.get_conduction_band_edge();
    double Efn = -dd.get_electron_electro_chemical_potential() + dd.get_electric_potential();
    double dE_n = 0;
    if (Efn <= Ec)
    {
      if (Efn >= Et)
        dE_n = Ec - Efn;
      else
        dE_n = Ec - Et;
    }
    double gamman = _tat->get_gamma(dd.get_electric_field().norm() * 100, T, dE_n);
    gamman = 1.0 / (gamman + 1);
    tau_n *= gamman;

    double Ev = dd.get_valence_band_edge();
    double Efp = -dd.get_hole_electro_chemical_potential() + dd.get_electric_potential();
    double dE_p = 0;
    if (Efp >= Ev)
    {
      if (Efp <= Et)
        dE_p = Efp - Ev;
      else
        dE_p = Et- Ev;
    }
    double gammap = _tat->get_gamma(dd.get_electric_field().norm() * 100, T, dE_p);
    gammap = 1.0 / (gammap + 1);
    tau_p *= gammap;
    //std::cerr << gamman << " "  << gammap << std::endl;
  }

  // TODO should take carrier temperatures
  double kT_e = T;
  double kT_h = T;
  double arg_e = get_trap_level() - Efn - dd.get_electric_potential();
  double arg_h = get_trap_level() - Efp - dd.get_electric_potential();
  double f_e, f_h;

  // get the occupations
  if (_dos == NULL)
  {
    std::pair<double, double> occ_e(Distributions::fermi_dirac(-arg_e, kT_e));
    f_e = occ_e.first;
    //deriv_e = occ_e.second;

    std::pair<double, double> occ_h(Distributions::fermi_dirac(-arg_h, kT_h));
    f_h = occ_h.first;
    //deriv_h = occ_h.second;
  }
  else
  {
    std::pair<double, double> occ_h(_dos->get_occupied_density_and_derivative(Efp, dd.get_electric_potential() - get_trap_level(), kT_h));
    f_h = occ_h.first;
    //f_h = _dos->get_occupied_density(-arg_h, kT_h);
    //deriv_h = _dos->get_occupied_density_derivative(-arg_h, kT_h);

    std::pair<double, double> occ_e(_dos->get_occupied_density_and_derivative(Efn, dd.get_electric_potential() - get_trap_level(), kT_e));
    f_e = occ_e.first;
    //f_e = _dos->get_occupied_density(-arg_e, kT_e);
    //deriv_e = _dos->get_occupied_density_derivative(-arg_e, kT_e);
  }

  double gc = (1.0 - f_e) / f_e;
  double gv = f_h / (1.0 - f_h);

  //if (arg_e / kT_e > 50)
    gc = exp(arg_e / kT_e);
  //if (arg_e / kT_e < -50)
    gv = exp(-arg_h / kT_h);



  double a = 1 + gc;
  double b = 1 + gv;
  double c = 1 - gc*gv;

  //a = 1 + exp((Et - Efn) / T);
  //b = 1 + exp((Efp - Et) / T);
  c = 1 - exp((Efp - Efn) / T);


  double denom = tau_p * n * a + tau_n * p * b + tau_n * tau_p * (_gen_VT + _gen_TC);
  double nom = n * p * c - (tau_p * n * gc * _gen_VT +
      tau_n * p * gv * _gen_TC + tau_n * tau_p * _gen_VT * _gen_TC);

  recomb_e = nom / denom;
  recomb_h = recomb_e;

}



void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double T = dd.get_lattice_temperature();

  double Et = get_trap_level();
  double f = std::exp((Et - dd.get_equilibrium_fermi_level()) / T);

  double dens = _density;
  if (_profile != nullptr)
    dens *= _profile->get_data(dd.get_element(), dd.get_coordinates());

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    dens = max(dens, 1e-32);
    tau_n = 1.0 / (dens * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (dens * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  if (_tat != NULL)
  {
    // definition of integration limit
    double Ec = dd.get_conduction_band_edge();
    double Efn = -dd.get_electron_electro_chemical_potential() + dd.get_electric_potential();
    double dE_n = 0;
    if (Efn <= Ec)
    {
      if (Efn >= Et)
        dE_n = Ec - Efn;
      else
        dE_n = Ec - Et;
    }
    double gamman = _tat->get_gamma(dd.get_electric_field().norm() * 100, T, dE_n);
    gamman = 1.0 / (gamman + 1);
    tau_n *= gamman;

    double Ev = dd.get_valence_band_edge();
    double Efp = -dd.get_hole_electro_chemical_potential() + dd.get_electric_potential();
    double dE_p = 0;
    if (Efp >= Ev)
    {
      if (Efp <= Et)
        dE_p = Efp - Ev;
      else
        dE_p = Et - Ev;
    }
    double gammap = _tat->get_gamma(dd.get_electric_field().norm() * 100, T, dE_p);
    gammap = 1.0 / (gammap + 1);
    tau_p *= gammap;
  }


  // TODO should take carrier temperatures
  double kT_e = T;
  double kT_h = T;
  double arg_e = get_trap_level() - Efn - dd.get_electric_potential();
  double arg_h = get_trap_level() - Efp - dd.get_electric_potential();
  double f_e, f_h;
  double deriv_e, deriv_h;

  //if (a < 0) a = p / denom;
  //if (b < 0) b = n / denom;
  // get the occupations
  if (_dos == NULL)
  {
    std::pair<double, double> occ_e(Distributions::fermi_dirac(-arg_e, kT_e));
    f_e = occ_e.first;
    deriv_e = occ_e.second;

    std::pair<double, double> occ_h(Distributions::fermi_dirac(-arg_h, kT_h));
    f_h = occ_h.first;
    deriv_h = occ_h.second;
  }
  else
  {
    std::pair<double, double> occ_h(_dos->get_occupied_density_and_derivative(Efp, dd.get_electric_potential() - get_trap_level(), kT_h));
    f_h = occ_h.first;
    deriv_h = occ_h.second;

    //f_h = _dos->get_occupied_density(-arg_h, kT_h);
    //deriv_h = _dos->get_occupied_density_derivative(-arg_h, kT_h);

    std::pair<double, double> occ_e(_dos->get_occupied_density_and_derivative(Efn, dd.get_electric_potential() - get_trap_level(), kT_e));
    f_e = occ_e.first;
    deriv_e = occ_e.second;
    //f_e = _dos->get_occupied_density(-arg_e, kT_e);
    //deriv_e = _dos->get_occupied_density_derivative(-arg_e, kT_e);
  }

  long double gc = (1.0 - f_e) / f_e;
  long double gv = f_h / (1.0 - f_h);
  long double deriv_gc = -deriv_e / (f_e * f_e);
  long double deriv_gv =  deriv_h / ((1 - f_h) * (1 - f_h));

  //if (arg_e / kT_e > 50)
  {
    gc = exp(arg_e / kT_e);
    deriv_gc = -gc / kT_e;
  }
  //if (arg_h / kT_h < -50)
  {
    gv = exp(-arg_h / kT_h);
    deriv_gv = gv / kT_h;
  }

  long double a = 1 + gc;
  long double b = 1 + gv;
  long double c = 1 - gc*gv;

  //a = 1 + exp((Et - Efn) / T);
  //b = 1 + exp((Efp - Et) / T);
  c = 1 - exp((Efp - Efn) / T);

  long double nom = n * p * c - (tau_p * n * gc * _gen_VT +
      tau_n * p * gv * _gen_TC + tau_n * tau_p * _gen_VT * _gen_TC);

  long double denom = tau_p * n * a + tau_n * p * b + tau_n * tau_p * (_gen_VT + _gen_TC);
  long double denom2 = denom * denom;

  long double r_e = nom / denom;
  long double r_h = r_e;

  long double dRedn = (p * c - tau_p * gc * _gen_VT - r_e * tau_p * a) / denom;
  long double dRedp = (n * c - tau_n * gv * _gen_TC - r_e * tau_n * b) / denom;


  long double dRedEfn = n * deriv_gc * (p * gv + tau_p * (_gen_VT + r_e)) / denom;
  long double dRedEfp = p * deriv_gv * (n * gc + tau_n * (_gen_TC + r_e)) / denom;


  recomb_e[0] = recomb_h[0] = dRedn;
  recomb_e[1] = recomb_h[1] = dRedp;
  recomb_e[2] = recomb_h[2] = dRedEfn;
  recomb_e[3] = recomb_h[3] = dRedEfp;
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

