// $Id$

#include "AugerBare.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "TiberMath.h"

#include "TiberModule.h"


void
AugerBare::read_database(void)
{
  /*
  const Database& db = get_database();
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

  */
}



void
AugerBare::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 3)
    throw InitFailedException("Bare Auger model needs exactly "
        "three recombining carriers");

  get_parameter("rate_constant", _rate_constant);
}




void
AugerBare::calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials)
{
  const ID id1 = this->get_carrier_ids()[0]; // electron
  const ID id2 = this->get_carrier_ids()[1]; // hole
  const ID id3 = this->get_carrier_ids()[2]; // final electron

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Efn = -dd.get_q_fermi_potential(id1);
  double Efp = -dd.get_q_fermi_potential(id2);
  double Efnf = -dd.get_q_fermi_potential(id3);
  double kT = dd.get_lattice_temperature();

  double n  = dd.get_q_density(id1);
  double p  = dd.get_q_density(id2);
  double dn  = dd.get_q_density_derivative(id1);
  double dp  = dd.get_q_density_derivative(id2);
  double qn = dd.get_carrier_properties(id1)->get_charge();
  double qp = dd.get_carrier_properties(id2)->get_charge();

  double Ecf = dd.get_carrier_properties(id3)->get_band_edge();
  auto   ff = Distributions::fermi_dirac(Efnf - Ecf, kT);

  double exponential = exp((Efp - 2*Efn + Efnf) / kT);
  double stat_fac = 1.0 - exponential;

  double g = _rate_constant * n * n * p * (1 - ff.first);

  double rec = g * stat_fac;

  R[id1] =  rec;
  R[id2] =  rec;
  R[id3] = -rec;

  // w.r.t. quasi Fermi levels
  double dR1 = _rate_constant * 2 * n * dn * p * (1 - ff.first) * stat_fac
      + 2 / kT * g * exponential;
  double dR2 = _rate_constant * n * n * dp * (1 - ff.first) * stat_fac
      - g / kT * exponential;
  double dR3 = _rate_constant * n * n * p * ff.second * stat_fac
      - g / kT * exponential;

  double dR0 = _rate_constant * stat_fac * (n * (2 * dn * p + n * dp) * (1 - ff.first)
      - n * n * p * ff.second);

  dPotentials[id1][id1] = -dR1;
  dPotentials[id1][id2] = -dR2;
  dPotentials[id1][id3] = -dR3;
  dPotentials[id2][id1] = -dR1;
  dPotentials[id2][id2] = -dR2;
  dPotentials[id2][id3] = -dR3;
  dPotentials[id3][id1] = dR1;
  dPotentials[id3][id2] = dR2;
  dPotentials[id3][id3] = dR3;
  dPotentials[id1][dd.n_known_carriers()] = dR0;
  dPotentials[id2][dd.n_known_carriers()] = dR0;
  dPotentials[id3][dd.n_known_carriers()] = -dR0;
}



