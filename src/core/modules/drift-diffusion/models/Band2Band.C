// $Id$

#include "Band2Band.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"

#include "TiberModule.h"


using namespace std;


void
Band2Band::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/bbt");

  _B_param = db.get("B", _B_param);
  _sigma = db.get("sigma", _sigma);
  _E0 = db.get("E0", _E0);

}



void
Band2Band::do_init(void)
{
  get_parameter("B", _B_param);
  get_parameter("sigma", _sigma);
  get_parameter("E0", _E0); // V/cm
}



void
Band2Band::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double gn = dd.get_electron_gamma();
  double gp = dd.get_hole_gamma();
  double E = dd.get_electric_field().size() + 1e-3;

  double D = (n * p - ni * ni * gn * gp) / ((p + ni) * (n + ni));
  recomb_e = recomb_h = _B_param * D * pow(E, _sigma) * exp(-_E0 / E);
}



void
Band2Band::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double gn = dd.get_electron_gamma();
  double gp = dd.get_hole_gamma();
  double E = dd.get_electric_field().size() + 1e-3;

  double D = (n * p - ni * ni * gn * gp) / ((p + ni) * (n + ni));
  double recomb = _B_param * pow(E, _sigma) * exp(-_E0 / E);

  recomb_e[0] = recomb_h[0] = recomb * ni / ((n + ni) * (n + ni)); // dR/dn
  recomb_e[1] = recomb_h[1] = recomb * ni / ((p + ni) * (p + ni)); // dR/dp
}


