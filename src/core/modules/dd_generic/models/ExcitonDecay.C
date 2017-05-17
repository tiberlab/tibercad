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


void
ExcitonDecay::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonDecay::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_ids().size() != 1)
    throw InitFailedException("ExcitonDecay model needs exactly one carrier");

  _tau = get_option("tau", _tau);

}




void
ExcitonDecay::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();


  double kT = dd.get_lattice_temperature();

  double id = get_carrier_ids()[0];
  unsigned int np = dd.n_known_carriers();

  double x = dd.get_q_density(id);
  double dx = dd.get_q_density_derivative(id);
  double Efx = -dd.get_q_fermi_potential(id);

  double beta = 1.0 / kT;
  double exponential = exp(-beta * Efx);
  double stat = 1.0 - exponential;

  double rate = stat * x / _tau;

  R[id] = rate;

  dPotentials[id][id] = - (beta * x * exponential + dx * stat) / _tau;

}
