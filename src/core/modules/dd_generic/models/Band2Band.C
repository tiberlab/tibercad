// $Id: Band2Band.C 3542 2013-03-01 09:31:59Z maufder $

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
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Band-to-band recombination model needs exactly "
        "two recombining carriers");

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  if (dd.get_carrier_properties(id1)->get_charge() *
      dd.get_carrier_properties(id2)->get_charge() == 0.0)
  {
    throw InitFailedException("Band-to-band recombination model implementation "
        "does not make sense for uncharged particles.");
  }

  get_parameter("B", _B_param);
  get_parameter("sigma", _sigma);
  get_parameter("E0", _E0); // V/cm
}



void
Band2Band::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double E01 = dd.get_carrier_properties(id1)->get_band_edge();
  double E02 = dd.get_carrier_properties(id2)->get_band_edge();
  double F = dd.get_electric_field().norm() + 1e-12;
  double dE = fabs(E01 - E02);

  double factor = _B_param * exp(-_E0 / F);
  //cerr << F << "  " << factor << endl;

  double q1 = dd.get_carrier_properties(id1)->get_charge();
  double q2 = dd.get_carrier_properties(id2)->get_charge();

  double kT = dd.get_lattice_temperature();

  if (q1 * q2 < 0)
  {

    if (q1 > 0)
    {
      swap(id1, id2);
      swap(q1, q2);
    }
    // now electron is id1

    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);
    double beta = 1.0/kT;

    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double dn1  = dd.get_q_density_derivative(id1);
    double dn2  = dd.get_q_density_derivative(id2);
    double q1 = dd.get_carrier_properties(id1)->get_charge();
    double q2 = dd.get_carrier_properties(id2)->get_charge();

    double exponential = exp((Ef2 - Ef1) * beta);
    double stat_fac = 1.0 - exponential;
    double g = factor * n1 * n2;

    R[id1] = g * stat_fac;
    R[id2] = g * stat_fac;
    //cerr << "n1 = " << n1 << " n2 = " << n2 << " -> " << R[id1] << endl;

    double dR0 = stat_fac * factor * (n2 * dn1 + n1 * dn2);
    double dR1 = -factor * n2 * (dn1 * stat_fac + beta * n1 * exponential);
    double dR2 = -factor * n1 * (dn2 * stat_fac - beta * n2 * exponential);

    dPotentials[id1][id1] = dR1;
    dPotentials[id1][id2] = dR2;
    dPotentials[id2][id1] = dR1;
    dPotentials[id2][id2] = dR2;
    dPotentials[id1][dd.n_known_carriers()] = dR0;
    dPotentials[id2][dd.n_known_carriers()] = dR0;
  }
  else
  {
    if (((q1 <= 0) && (E01 > E02)) ||
        ((q1 > 0)  && (E01 < E02)))
    {
      swap(id1, id2);
    }

    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double N1  = dd.get_carrier_properties(id1)->get_maximum_density();
    double N2  = dd.get_carrier_properties(id2)->get_maximum_density();
    double dn1 = dd.get_q_density_derivative(id1);
    double dn2 = dd.get_q_density_derivative(id2);
    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);
    double beta = -dd.get_carrier_properties(id1)->get_charge_sign();

    double exponential1 = exp((Ef2 - Ef1) * beta);
    double exponential2 = exp((Ef1 - Ef2) * beta);
    double stat1 = 1.0 - exponential1;
    double stat2 = 1.0 - exponential2;

    double n1f = n1/N1;
    double n2f = n2/N2;

    // TODO N2 should be the maximum density for n2, but now its 100*Nc
    R[id1] = factor * (stat1 * n1 * (1 - n2f) - stat2 * n2 * (1 - n1f));
    R[id2] = -R[id1];


    double dR0 =  factor * (stat1 * ( (1 - n2f)*dn1 - n1*dn2/N2 )
        - stat2 * ( (1 - n1f)*dn2 - n2 * dn1/N1 ));
    double dR1 = -factor * ((1 - n2f) * ( dn1 * stat1 + beta * n1 * exponential1)
        - n2 * ( dn1/N1 * stat2 + beta * (1 - n1f) * exponential2));
    double dR2 = -factor * (n1 * (-dn2/N2 * stat1 - beta * (1 - n2f) * exponential1)
        + (1 - n1f) * (dn2 * stat2 + beta * n2 * exponential2));

    dPotentials[id1][id1] =  dR1;
    dPotentials[id1][id2] =  dR2;
    dPotentials[id2][id1] = -dR1;
    dPotentials[id2][id2] = -dR2;
    dPotentials[id1][dd.n_known_carriers()] =  dR0;
    dPotentials[id2][dd.n_known_carriers()] = -dR0;
  }
}
