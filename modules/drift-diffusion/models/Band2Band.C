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

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double E = dd.get_electric_field().norm() + 1e-3;
  double T = dd.get_lattice_temperature();

  double expf = exp((Efp - Efn) / T);
  double c = 1.0 - expf;

  double D = n * p / ((p + ni) * (n + ni));
  recomb_e = recomb_h = _B_param * c * D * pow(E, _sigma) * exp(-_E0 / E);
}



void
Band2Band::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double E = dd.get_electric_field().norm() + 1e-3;
  double T = dd.get_lattice_temperature();

  double expf = exp((Efp - Efn) / T);
  double c = 1.0 - expf;

  double D = n * p / ((p + ni) * (n + ni));
  double factor = _B_param * pow(E, _sigma) * exp(-_E0 / E);
  double recomb = c * D * factor;

  double tmp = p / ((p + ni) * (n + ni)) - D / (n + ni);
  recomb_e[0] = recomb_h[0] = tmp * c * factor; // dR/dn
  tmp = n / ((p + ni) * (n + ni)) - D / (p + ni);
  recomb_e[1] = recomb_h[1] = tmp * c * factor; // dR/dp
  recomb_e[2] = recomb_h[2] = -D * factor / T * expf;
  recomb_e[3] = recomb_h[3] = D * factor / T * expf;
}


