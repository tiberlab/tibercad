// $Id: ExcitonGeneration.C 3542 2013-03-01 09:31:59Z maufder $

#include "ExcitonGeneration.h"
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
ExcitonGeneration::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonGeneration::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_ids().size() != 3)
  {
    throw InitFailedException("ExcitonGeneration model requires three carriers in input.");
  }

  //vector<string> carriers;
  //get_option("carriers", carriers);
  //reorder_ids(carriers);
  _gamma = get_option("gamma", _gamma);
  _gamma = get_option("C", _gamma);


  _stat_fac = get_option("stat_fac", _stat_fac);

}




void
ExcitonGeneration::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double beta = 1.0 / dd.get_lattice_temperature();


  // for the moment we rely on the ordering of the carriers!

    ID idn = get_carrier_ids()[0];
    ID idp = get_carrier_ids()[1];
    ID idx = get_carrier_ids()[2];

    char ct1 = dd.get_carrier_properties(idn)->get_carrier_type();
    char ct2 = dd.get_carrier_properties(idp)->get_carrier_type();

    if (ct1 == 'h')
      swap(idn, idp);

    unsigned int nc = dd.n_known_carriers();

    double spin = dd.get_carrier_properties(idx)->get_spin();
    double fac = (spin == 0.0) ? 0.25 : 0.75;

    if (!_stat_fac)
      fac = 1.0;

    double x = dd.get_q_density(idx);
    double n = dd.get_q_density(idn);
    double p = dd.get_q_density(idp);
    double dx = dd.get_q_density_derivative(idx);
    double dn = dd.get_q_density_derivative(idn);
    double dp = dd.get_q_density_derivative(idp);
    double Nx = dd.get_carrier_properties(idx)->get_effective_DOS();
    double Efx = -dd.get_q_fermi_potential(idx);
    double Efn = -dd.get_q_fermi_potential(idn);
    double Efp = -dd.get_q_fermi_potential(idp);

    double exponential = exp(beta*(Efx - Efn + Efp));
    double stat = 1.0 - exponential;

    double rate = _gamma * stat * n * p * ((Nx + x) / Nx);

    R[idx] = -rate;
    R[idn] = rate;
    R[idp] = rate;

    double derf = _gamma * stat * ((Nx + x) / Nx) * (n * dp + p * dn);
    double derx = _gamma * n * p * (beta * exponential * ((Nx + x) / Nx) + (dx / Nx) * stat);
    double dern = -_gamma * p * (dn * stat + beta * n * exponential) * ((Nx + x) / Nx);
    double derp = -_gamma * n * (dp * stat - beta * p * exponential) * ((Nx + x) / Nx);

    dPotentials[idx][idx] = - derx;
    dPotentials[idx][idn] = - dern;
    dPotentials[idx][idp] = - derp;
    dPotentials[idx][nc]  = - derf;

    dPotentials[idn][idx] = derx;
    dPotentials[idn][idn] = dern;
    dPotentials[idn][idp] = derp;
    dPotentials[idn][nc]  = derf;

    dPotentials[idp][idx] = derx;
    dPotentials[idp][idn] = dern;
    dPotentials[idp][idp] = derp;
    dPotentials[idp][nc]  = derf;

}


