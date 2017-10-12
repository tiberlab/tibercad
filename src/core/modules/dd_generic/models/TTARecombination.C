// $Id$

#include "TTARecombination.h"
#include "DriftDiffusionProperties.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;

namespace {
  pair<double, double> fermi_dirac(double E)
  {
    double f = 0, deriv = 0;
    double arg = E;
    if (arg > 50)
    {
      f = exp(-arg);
      deriv = -f;
    }
    else if (arg < -50)
    {
      deriv = -exp(arg);
      f = 1 + deriv;
    }
    else
    {
      double expfac = exp(arg);
      double denom = 1.0 + expfac;
      f = 1.0 / denom;
      deriv = -expfac * f / denom;
    }

    return make_pair(f, deriv);
  }
}

/*
void
TTARecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);
}
*/



void
TTARecombination::do_init(void)
{
  RecombinationModelInterface::do_init();

  get_parameter("C", C_);

  if (get_carrier_names().size() != 1)
    throw InitFailedException("Triplet-Triplet annihilation model "
        "needs specification of exactly one carrier");

  this->set_exponents(vector<double>(1, 2.0));

}




void
TTARecombination::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

  ID id1 = this->get_carrier_ids()[0];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();

  double kT = dd.get_lattice_temperature();

  double Ef1 = -dd.get_q_fermi_potential(id1);
  double beta = 1.0/kT;

  double n1  = dd.get_q_density(id1);
  double dn1  = dd.get_q_density_derivative(id1);
  double q1 = dd.get_carrier_properties(id1)->get_charge();


  double exponential = exp(-2 * Ef1 * beta);
  double stat_fac = 1.0 - exponential;
  double g = C_ * n1 * n1;

  R[id1] = g * stat_fac;

  double dR0 = stat_fac * C_ * 2 * n1 * dn1;
  double dR1 = -g * 2 * beta * exponential;

  dPotentials[id1][id1] = dR0 + dR1;
}
