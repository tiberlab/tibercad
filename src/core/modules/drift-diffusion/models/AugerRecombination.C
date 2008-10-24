// $Id$

#include "AugerRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"




TIBER_MODULE(AugerRecombination, auger)



void
AugerRecombination::read_database(void)
{
  Database& db = get_database();
  db.set_section("recombination/auger");

  std::vector<double> data(2);
  if (db.has_variable("A"))
  {
    db.get("A", data);
    if (data.size() < 2)
      throw InitFailedException("Auger recombination needs two "
          "values for \'A\' in database.");
    _An = data[0];
    _Ap = data[1];
  }

  if (db.has_variable("B"))
  {
    db.get("B", data);
    if (data.size() < 2)
      throw InitFailedException("Auger recombination needs two "
          "values for \'B\' in database.");
    _Bn = data[0];
    _Bp = data[1];
  }

  if (db.has_variable("C"))
  {
    db.get("C", data);
    if (data.size() < 2)
      throw InitFailedException("Auger recombination needs two "
          "values for \'C\' in database.");
    _Cn = data[0];
    _Cp = data[1];
  }

  if (db.has_variable("H"))
  {
    db.get("H", data);
    if (data.size() < 2)
      throw InitFailedException("Auger recombination needs two "
          "values for \'H\' in database.");
    _Hn = data[0];
    _Hp = data[1];
  }

  if (db.has_variable("N0"))
  {
    db.get("N0", data);
    if (data.size() < 2)
      throw InitFailedException("Auger recombination needs two "
          "values for \'N0\' in database.");
    _N0n = data[0];
    _N0p = data[1];
  }

}



void
AugerRecombination::do_init(void)
{

  if (has_parameter("Cn"))
  {
    _Cn = get_parameter("Cn", _Cn);
    _fixed_Cn = true;
  }
  if (has_parameter("Cp"))
  {
    _Cp = get_parameter("Cp", _Cp);
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
  long double ni2 = ni * ni;
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
  long double ni2 = ni * ni;
  long double np = n * p;

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * (np - ni2);
  long double B = Cp * (np - ni2);

  recomb_e[0] = recomb_h[0] = Cn * np + Cp * p * p + A;
  recomb_e[1] = recomb_h[1] = Cp * np + Cn * n * n + B;
}



void
AugerRecombination::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const AugerRecombination* scA =
    dynamic_cast<const AugerRecombination*>(comp_A);
  const AugerRecombination* scB =
    dynamic_cast<const AugerRecombination*>(comp_B);

  _Cn = alloy(scA->_Cn, scB->_Cn, xa);
  _Cp = alloy(scA->_Cp, scB->_Cp, xa);
}

