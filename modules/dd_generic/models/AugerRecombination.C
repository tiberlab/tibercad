// $Id: AugerRecombination.C 3880 2014-07-02 14:47:37Z maufder $

#include "AugerRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"

#include "TiberModule.h"


void
AugerRecombination::read_database(void)
{
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

}



void
AugerRecombination::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Auger recombination model needs exactly "
        "two recombining carriers");

  // look for electron and hole

  if (has_parameter("Cn"))
  {
    get_parameter("Cn", _Cn);
    _fixed_Cn = true;
  }
  if (has_parameter("Cp"))
  {
    get_parameter("Cp", _Cp);
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
    Cn = tmp * std::exp(-dd.get_q_density(this->get_carrier_ids()[0]) / _N0n);
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
    Cp = tmp * std::exp(-dd.get_q_density(this->get_carrier_ids()[1]) / _N0p);
  }
  return Cp;
}


void
AugerRecombination::calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials)
{
  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  //double E01 = dd.get_carrier_properties(id1)->get_band_edge();
  //double E02 = dd.get_carrier_properties(id2)->get_band_edge();

  const char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();
  const char ct2 = dd.get_carrier_properties(id2)->get_carrier_type();

  if (ct1 != ct2)
  {
    // this is the standard direct recombination for electrons hole pairs
    if (ct1 != 'e')
      std::swap(id1, id2);

    double Efn = -dd.get_q_fermi_potential(id1);
    double Efp = -dd.get_q_fermi_potential(id2);
    double kT = dd.get_lattice_temperature();

    double n  = dd.get_q_density(id1);
    double p  = dd.get_q_density(id2);
    double dn  = dd.get_q_density_derivative(id1);
    double dp  = dd.get_q_density_derivative(id2);
    double qn = dd.get_carrier_properties(id1)->get_charge();
    double qp = dd.get_carrier_properties(id2)->get_charge();

    double exponential = exp((Efp - Efn) / kT);
    double stat_fac = 1.0 - exponential;

    double Cn = get_Cn();
    double Cp = get_Cp();

    double A = Cn * n;
    double B = Cp * p;
    double Rec = (A + B) * n * p * stat_fac;

    R[id1] = R[id2] = Rec;

    double dRedn = (2*A + B) * p * stat_fac;
    double dRedp = (2*B + A) * n * stat_fac;
    double dRedEfn = -(A + B) * n * p * exponential / kT;
    double dRedEfp = -dRedEfn;
    double dR0 = -dRedn - dRedp;
    double dR1 = dRedn * dn + dRedEfn;
    double dR2 = dRedp * dp + dRedEfp;

    dPotentials[id1][id1] = dR1;
    dPotentials[id1][id2] = dR2;
    dPotentials[id2][id1] = dR1;
    dPotentials[id2][id2] = dR2;
    dPotentials[id1][dd.n_known_carriers()] = dR0;
    dPotentials[id2][dd.n_known_carriers()] = dR0;
  }

}

