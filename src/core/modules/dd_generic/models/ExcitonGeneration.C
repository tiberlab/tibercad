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


  _stat_fac = get_option("stat_fac", _stat_fac);


  if (_C.size() == 1)
    _C.resize(carriers.size(), _C[0]);

  if ( (_C.size() > 1) && (_C.size() != carriers.size()) )
    throw InitFailedException("Number of excitons not consistent with number of C's");

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  for (auto name : get_carrier_names())
  {
    if (!dd.get_carrier_properties(name)->is_exciton())
      throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");

    double spin = dd.get_carrier_properties(name)->get_spin();
    if ((spin != 0.0) && (spin != 1.0))
      throw InitFailedException("Recombination '" + get_default_name() + ": '" + name + "' spin can be either 0 or 1");
  }
  */

  /*
  for (auto id : get_carrier_ids())
  {
    _exciton_carriers.insert( make_pair(id, vector<unsigned int>()) );
    vector<string> ex_carriers = dd.get_carrier_properties(id)->get_exciton_carriers();
    for (auto exc : ex_carriers)
      _exciton_carriers[id].push_back( dd.get_carrier_id(exc) );
  }
  */



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

    double rate = _C[i] * stat * n * p * (Nx + x) / Nx;

    R[idx] = -rate;
    R[idn] = rate;
    R[idp] = rate;

    double derf = _C[i] * stat * (Nx + x) * (n * dp + p * dn) / Nx;
    double derx = - _C[i] * n * p * (dx * stat - beta * (Nx + x) * exponential) / Nx;
    double dern = - _C[i] * p * (Nx + x) * (dn * stat + beta * n * exponential) / Nx;
    double derp = - _C[i] * n * (Nx + x) * (dp * stat - beta * n * exponential) / Nx;

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


