// $Id$

#include "AugerRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"




TIBER_MODULE(AugerRecombination, recombination, auger)



void
AugerRecombination::read_database(void)
{
  Database& db = get_database();
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
    Cn = tmp * std::exp(-dd.get_electron_density() / _N0n);
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
    Cp = tmp * std::exp(-dd.get_electron_density() / _N0p);
  }
  return Cp;
}



void
AugerRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double ni = dd.get_intrinsic_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  long double ni2 = ni * ni * gn * gp;
  long double np = n * p;

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * n;
  long double B = Cp * p;
  recomb_e = recomb_h = (A + B) * (np - ni2);
}



void
AugerRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double ni = dd.get_intrinsic_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  long double ni2 = ni * ni * gn * gp;
  long double np = n * p;

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * (np - ni2);
  long double B = Cp * (np - ni2);

  recomb_e[0] = recomb_h[0] = Cn * np + Cp * p * p + A;
  recomb_e[1] = recomb_h[1] = Cp * np + Cn * n * n + B;
}




