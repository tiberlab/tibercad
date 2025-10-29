// $Id: EPQuenching.C 4492 2017-11-23 13:42:01Z maufder $

#include "EPQuenching.h"
#include "DriftDiffusionProperties.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;


/*
void
EPQuenching::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);
}
*/



void
EPQuenching::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 1)
    throw InitFailedException("Triplet-Polaron quenching model "
        "needs specification of exactly one carrier (triplet)");

  get_parameter("C", C_);
  string quencher;
  get_parameter("quencher", quencher);

  _quencher = get_driftdiffusionproperties().get_carrier_id(quencher);
  if (_quencher < 0)
  {
    throw InitFailedException("Triplet-Polaron quenching model "
        "needs specification of the quenching particle 'quencher'");
  }


}




void
EPQuenching::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

  ID id1 = this->get_carrier_ids()[0];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();

  double kT = dd.get_lattice_temperature();

  double Ef1 = -dd.get_q_fermi_potential(id1);
  double beta = 1.0/kT;

  double n1  = dd.get_q_density(id1);
  double dn1 = dd.get_q_density_derivative(id1);
  double q1  = dd.get_carrier_properties(id1)->get_charge();

  double nq  = dd.get_q_density(_quencher);
  double dnq = dd.get_q_density_derivative(_quencher);
  double qq  = dd.get_carrier_properties(_quencher)->get_charge();

  double exponential = exp(-Ef1 * beta);
  double stat_fac = 1.0 - exponential;
  double g = C_ * n1 * nq;

  R[id1] = g * stat_fac;

  double dR0 = -stat_fac * C_ * nq * dn1;
  double dR1 = -g *  beta * exponential;
  double dRq = -stat_fac * C_ * (n1 * dnq);

  dPotentials[id1][id1] = dR0 + dR1;
  dPotentials[id1][_quencher] = dRq;
  dPotentials[id1][dd.n_known_carriers()] = -fabs(qq)*dRq;
}
