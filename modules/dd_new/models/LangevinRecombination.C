

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
}



double
LangevinRecombination::calculate_rate_and_derivatives(std::vector<double>& dPotentials)
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

  double mun = dd.get_q_mobility(id1);
  double mup = dd.get_q_mobility(id2);

  double exponential = exp((Efp - Efn) / kT);
  double stat_fac = 1.0 - exponential;
  double prefactor = _gamma * Constants::e * 100 / (_er * Constants::e0) * (mun + mup);

  double g = prefactor * n * p;

  double R = g * stat_fac;

  double dR1 = prefactor * p * (dn * stat_fac - 1/kT * n * exponential);
  double dR2 = prefactor * n * (dp * stat_fac + 1/kT * p * exponential);

  dPotentials[id1] = dR1;
  dPotentials[id2] = dR2;
  dPotentials[dd.n_known_carriers()] = stat_fac * prefactor * (p * dn + n * dp);

  return R;
}

