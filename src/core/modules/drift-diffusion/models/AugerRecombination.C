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

  _Cn = db.get("Cn", _Cn);
  _Cp = db.get("Cp", _Cp);
}



void
AugerRecombination::do_init(void)
{
  _Cn = get_parameter("Cn", _Cn);
  _Cp = get_parameter("Cp", _Cp);
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

  long double A = _Cn * n * (np - ni2);
  long double B = _Cp * p * (np - ni2);
  recomb_e = recomb_h = A + B;
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

  long double A = _Cn * (np - ni2);
  long double B = _Cp * (np - ni2);

  recomb_e[0] = recomb_h[0] = _Cn * np + _Cp * p * p + A;
  recomb_e[1] = recomb_h[1] = _Cp * np + _Cn * n * n + B;
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

