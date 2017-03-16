// $Id: ExcitonISC.C 3542 2013-03-01 09:31:59Z maufder $

#include "ExcitonISC.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;

ExcitonISC::QRecMap
ExcitonISC::_qrec_vals;

void
ExcitonISC::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonISC::do_init(void)
{
  RecombinationModelInterface::do_init();


  if (get_carrier_names().size() != 2)
    throw InitFailedException("Exciton ISC model needs exactly "
        "two recombining carriers");

  get_parameter("C", _C);

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  vector<double> spin;
  for (auto name : get_carrier_names())
  {
    if (!dd.get_carrier_properties(name)->is_exciton())
      throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");

    spin.push_back( dd.get_carrier_properties(name)->get_spin() );
  }

  if (!( ((spin[0] == 1.0)&&(spin[1] == 0.0)) || ((spin[0] == 0.0)&&(spin[1] == 1.0)) ))
    throw InitFailedException("Recombination '" + get_default_name() + ": excitons must be a singlet and a triplet");


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
ExcitonISC::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  double E1  = dd.get_carrier_properties(id1)->get_band_edge();
  double E2  = dd.get_carrier_properties(id2)->get_band_edge();

  if (E1<E2)
  {
    swap(id1, id2);
    swap(E1, E2);
  }

  double kT = dd.get_lattice_temperature();
  double x1  = dd.get_q_density(id1);
  double x2  = dd.get_q_density(id2);
  double N2  = dd.get_carrier_properties(id2)->get_effective_DOS();
  double dx1 = dd.get_q_density_derivative(id1);
  double dx2 = dd.get_q_density_derivative(id2);
  double f1 = dd.get_q_fermi_potential(id1);
  double f2 = dd.get_q_fermi_potential(id2);

  double beta = 1.0/kT;
  double exponential = exp( beta*(f1-f2) );
  double stat = 1.0 - exponential;

  double rate = _C * stat * x1 * (N2 + x2);
  double der1 = - _C * (N2 + x2) * (dx1 * stat + beta * x1 * exponential);
  double der2 = - _C * x1 * (dx2 * stat + beta * (N2 + x2) * exponential);

  R[id1] = rate;
  R[id2] = -rate;

  dPotentials[id1][id1] =  der1;
  dPotentials[id1][id2] =  der2;
  dPotentials[id2][id1] = -der1;
  dPotentials[id2][id2] = -der2;

}

void
ExcitonISC::do_reinit(void)
{

}
