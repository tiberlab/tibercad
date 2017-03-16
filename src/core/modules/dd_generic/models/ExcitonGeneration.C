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

ExcitonGeneration::QRecMap
ExcitonGeneration::_qrec_vals;

void
ExcitonGeneration::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonGeneration::do_init(void)
{
  RecombinationModelInterface::do_init();


  vector<string> carriers;

  get_option("carriers", carriers);
  reorder_ids(carriers);
  get_option("C", _C);


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

  for (auto id : get_carrier_ids())
  {
    _exciton_carriers.insert( make_pair(id, vector<unsigned int>()) );
    vector<string> ex_carriers = dd.get_carrier_properties(id)->get_exciton_carriers();
    for (auto exc : ex_carriers)
      _exciton_carriers[id].push_back( dd.get_carrier_id(exc) );
  }


  string quantumsim = get_option("optics_simulation", "");
  if (!quantumsim.empty())
  {
    _quantum_optics = SimulationInterface::find_simulation(quantumsim);
    if (_quantum_optics == NULL)
      throw InitFailedException("Cannot find optics simulation \'" + quantumsim + "\'");

    _rec_id = _quantum_optics->get_solution_id("Recombination");
    if (_rec_id == INVALID_ID)
      throw InitFailedException("Simulation \'" + quantumsim + "\'" +
          " does not have the needed solution \'Recombination\'");

    // TODO should check for consistency of regions
  }

}




void
ExcitonGeneration::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double beta = 1.0 / dd.get_lattice_temperature();

  for (size_t i = 0; i<_C.size(); i++)
  {
    ID idx = get_carrier_ids()[i];
    ID id1 = _exciton_carriers[idx][0];
    ID id2 = _exciton_carriers[idx][1];
    char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();
    char ct2 = dd.get_carrier_properties(id2)->get_carrier_type();

    ID idn = (ct1 == 'e') ? id1 : id2;
    ID idp = (ct2 == 'h') ? id2 : id1;

    unsigned int nc = dd.n_known_carriers();

    double spin = dd.get_carrier_properties(idx)->get_spin();
    double fac = (spin == 0.0) ? 0.25 : 0.75;

    double x = dd.get_q_density(idx);
    double n = dd.get_q_density(idn);
    double p = dd.get_q_density(idp);
    double dx = dd.get_q_density_derivative(idx);
    double dn = dd.get_q_density_derivative(idn);
    double dp = dd.get_q_density_derivative(idp);
    double Nx = dd.get_carrier_properties(idx)->get_effective_DOS();
    double fx = dd.get_q_fermi_potential(idx);
    double fn = dd.get_q_fermi_potential(idn);
    double fp = dd.get_q_fermi_potential(idp);

    double exponential = exp(beta*(fn-fp-fx));
    double stat = 1.0 - exponential;

    double rate = _C[i] * stat * n * p * (Nx + x);

    R[idx] = - fac * rate;
    R[idn] = rate;
    R[idp] = rate;

    double derf = _C[i] * stat * (Nx + x) * (n * dp + p * dn);
    double derx = - _C[i] * n * p * (dx * stat - beta * (Nx + x) * exponential);
    double dern = - _C[i] * p * (Nx + x) * (dn * stat + beta * n * exponential);
    double derp = - _C[i] * n * (Nx + x) * (dp * stat - beta * n * exponential);

    dPotentials[idx][idx] = - fac * derx;
    dPotentials[idx][idn] = - fac * dern;
    dPotentials[idx][idp] = - fac * derp;
    dPotentials[idx][nc]  = - fac * derf;

    dPotentials[idn][idx] = derx;
    dPotentials[idn][idn] = dern;
    dPotentials[idn][idp] = derp;
    dPotentials[idn][nc]  = derf;

    dPotentials[idp][idx] = derx;
    dPotentials[idp][idn] = dern;
    dPotentials[idp][idp] = derp;
    dPotentials[idp][nc]  = derf;
  }

}

void
ExcitonGeneration::do_reinit(void)
{

}
