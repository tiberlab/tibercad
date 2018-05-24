// $Id: AugerRecombination.C 3880 2014-07-02 14:47:37Z maufder $

#include "AugerRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"

#include "TiberModule.h"


void
AugerRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/auger");

  std::vector<double> data(2, 0);
  db.get("A", data, true);
  _An = data[0];
  _Ap = data[1];

  data = std::vector<double>(2, 0);
  db.get("B", data);
  _Bn = data[0];
  _Bp = data[1];

  data = std::vector<double>(2, 0);
  db.get("C", data);
  _Cn = data[0];
  _Cp = data[1];

  data = std::vector<double>(2, 0);
  db.get("H", data);
  _Hn = data[0];
  _Hp = data[1];

  data = std::vector<double>(2, 0);
  db.get("N0", data);
  _N0n = data[0];
  _N0p = data[1];

}



void
AugerRecombination::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Auger recombination model needs exactly "
        "two recombining carriers");

  // look for electron and hole

  if (has_parameter("Cn"))
  {
    get_parameter("Cn", _Cn);
    _fixed_Cn = true;
  }
  if (has_parameter("Cp"))
  {
    get_parameter("Cp", _Cp);
    _fixed_Cp = true;
  }
}


inline
double
AugerRecombination::get_Cn(void)
{
  double Cn = _Cn;
  if (!_fixed_Cn)
  {
    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
    double T = dd.get_lattice_temperature() / T0;
    double tmp = _An + _Bn * T + _Cn * T * T;
    Cn = tmp * std::exp(-dd.get_q_density(this->get_carrier_ids()[0]) / _N0n);
  }
  return Cn;
}

inline
double
AugerRecombination::get_Cp(void)
{
  double Cp = _Cp;
  if (!_fixed_Cp)
  {
    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
    double T = dd.get_lattice_temperature() / T0;
    double tmp = _Ap + _Bp * T + _Cp * T * T;
    Cp = tmp * std::exp(-dd.get_q_density(this->get_carrier_ids()[1]) / _N0p);
  }
  return Cp;
}


double
AugerRecombination::calculate_rate_and_derivatives(std::vector<double>& dPotentials)
{
  const ID id1 = this->get_carrier_ids()[0];
  const ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Efn = -dd.get_q_fermi_potential(id1);
  double Efp = -dd.get_q_fermi_potential(id2);
  double kT = dd.get_lattice_temperature();

  double n  = dd.get_q_density(id1);
  double p  = dd.get_q_density(id2);
  double dn  = dd.get_q_density_derivative(id1);
  double dp  = dd.get_q_density_derivative(id2);
  double qn = dd.get_carrier_properties(id1)->get_charge();
  double qp = dd.get_carrier_properties(id2)->get_charge();

  //double E01 = dd.get_carrier_properties(id1)->get_band_edge();
  //double E02 = dd.get_carrier_properties(id2)->get_band_edge();

  double exponential = exp((Efp - Efn) / kT);
  double stat_fac = 1.0 - exponential;
  /*
  double g = C_ * n1 * n2;

  double R = g * stat_fac;

  double dR1 = -C_ * n2 * (-dn1 * stat_fac + 1/kT * n1 * exponential);
  double dR2 = -C_ * n1 * (-dn2 * stat_fac - 1/kT * n2 * exponential);

  dPotentials[id1] = dR1;
  dPotentials[id2] = dR2;
  dPotentials[dd.n_known_carriers()] = stat_fac * C_ * (n2 * dn1 + n1 * dn2);
  */
  double R = 0;
  return R;
}

/*
void
AugerRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature();
  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double np = n * p;

  long double g = 1.0 - exp((Efp - Efn) / T);

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * n;
  long double B = Cp * p;
  recomb_e = recomb_h = (A + B) * np * g;
}



void
AugerRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature();
  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double np = n * p;

  long double e = exp((Efp - Efn) / T);
  long double g = 1 - e;

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * n;
  long double B = Cp * p;

  long double dRedn = (2*A + B) * p * g;
  long double dRedp = (2*B + A) * n * g;
  long double dRedEfn = -(A + B) * np * e / T;
  long double dRedEfp = -dRedEfn;

  recomb_e[0] = recomb_h[0] = dRedn;
  recomb_e[1] = recomb_h[1] = dRedp;
  recomb_e[2] = recomb_h[2] = dRedEfn;
  recomb_e[3] = recomb_h[3] = dRedEfp;
}
*/



