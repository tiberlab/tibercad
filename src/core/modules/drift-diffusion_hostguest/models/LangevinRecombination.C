

#include "LangevinRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"



using namespace std;


void
LangevinRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("permittivity");

  _er = db.get("permittivity", _er);
}

void
LangevinRecombination::do_init(void)
{
  get_parameter("gamma", _gamma);
  get_parameter("mindens", _mindens);
}

void
LangevinRecombination::get_net_recombination_rates(double& recomb_e, double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double mun = dd.get_electron_mobility();
  double mup = dd.get_hole_mobility();
  
  if (((n < _mindens) || (p < _mindens)) && ((n0 < _mindens) || (p0 < _mindens))) {
    recomb_e = recomb_h = 0.0;
  }
  else {
    recomb_e = recomb_h = _gamma * Constants::e * 100.0 / (_er * Constants::e0) * (mun + mup) * (n*p - n0*p0);
  }
}


void
LangevinRecombination::get_net_recombination_rate_derivatives(std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double mun = dd.get_electron_mobility();
  double mup = dd.get_hole_mobility();
  double dn_dphi = dd.get_electron_density_derivative();
  double dp_dphi = dd.get_hole_density_derivative();
  double dmun_dphi = dd.get_electron_mobility_derivative_potential();
  double dmup_dphi = dd.get_hole_mobility_derivative_potential();

  if (((n < _mindens) || (p < _mindens)) && ((n0 < _mindens) || (p0 < _mindens))) {
    recomb_e[0] = recomb_h[0] = 0.0;
    recomb_e[1] = recomb_h[1] = 0.0;
  }
  else {
    recomb_e[0] = recomb_h[0] = _gamma * Constants::e * 100.0 / (_er * Constants::e0) * ( (dmun_dphi / dn_dphi + dmup_dphi / dn_dphi) * (n*p - n0*p0) + (mun + mup) * p ); // dR/dn
    recomb_e[1] = recomb_h[1] = _gamma * Constants::e * 100.0 / (_er * Constants::e0) * ( (dmun_dphi / dp_dphi + dmup_dphi / dp_dphi) * (n*p - n0*p0) + (mun + mup) * n ); // dR/dp
  }
}

void
LangevinRecombination::do_reinit(void)
{
  
}
