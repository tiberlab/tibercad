// $Id: ExcitonDecay.C 3542 2013-03-01 09:31:59Z maufder $

#include "ExcitonDecay.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;

ExcitonDecay::QRecMap
ExcitonDecay::_qrec_vals;

void
ExcitonDecay::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonDecay::do_init(void)
{
  RecombinationModelInterface::do_init();


  vector<string> carriers;

  get_option("carriers", carriers);
  reorder_ids(carriers);
  get_option("tau", _tau);


  if (_tau.size() == 1)
    _tau.resize(carriers.size(), _tau[0]);

  if ( (_tau.size() > 1) && (_tau.size() != carriers.size()) )
    throw InitFailedException("Number of excitons not consistent with number of tau's");

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  for (auto name : get_carrier_names())
  {
    if (!dd.get_carrier_properties(name)->is_exciton())
      throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");

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
ExcitonDecay::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();


  double kT = dd.get_lattice_temperature();

  for (size_t i = 0; i<_tau.size(); i++)
  {
    double id = get_carrier_ids()[i];
    unsigned int np = dd.n_known_carriers();

    double x = dd.get_q_density(id);
    double dx = dd.get_q_density_derivative(id);
    double fx = dd.get_q_fermi_potential(id);

    double beta = 1.0 /kT;
    double exponential = exp(beta * fx);
    double stat = 1.0 - exponential;

    double rate = stat * x / _tau[i];

    R[id] = rate;

    dPotentials[id][id] = - (beta * x * exponential + dx * stat) / _tau[i];
  }

}

void
ExcitonDecay::do_reinit(void)
{

}
