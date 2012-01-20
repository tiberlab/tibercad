// $Id$

#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"




TIBER_MODULE(DirectRecombination, recombination, direct)

using namespace std;


void
DirectRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);

}



void
DirectRecombination::do_init(void)
{
  get_parameter("C", C_);

  string quantumsim = get_option("optics_simulation", "");
  if (!quantumsim.empty())
  {
    _quantum_optics = SimulationInterface::find_simulation(quantumsim);
    if (_quantum_optics == NULL)
      throw InitFailedException("Cannot find optics simulation \'" + quantumsim + "\'");

    _rec_id = _quantum_optics->get_solution_id("recombination");
    if (_rec_id == INVALID_ID)
      throw InitFailedException("Simulation \'" + quantumsim + "\'" +
          " does not have the needed solution \'recombination\'");
  }

}



void
DirectRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;

  recomb_e = recomb_h = C_ * (n * p - ni * ni * gn * gp);
}



void
DirectRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = C_ * p; // dR/dn
  recomb_e[1] = recomb_h[1] = C_ * n; // dR/dp
}


void
DirectRecombination::do_reinit(void)
{
  if (_quantum_optics != NULL)
  {
    map<ID, vector<double> > data;
    data[_rec_id];

    if (_quantum_optics->get_solution(data))
    {
      double rec = data[_rec_id][0];

      // now we have to integrate the term (np - ni^2)
      // for that, we have to loop over all elements
      // TODO this has to be checked for parallel execution
      const SimulationEnvironment& env =
          SimulationInterface::get_simulation(get_simulator_id())->get_environment();
      SimulationEnvironment::ConstElemIterator it(env.elements_begin());
      SimulationEnvironment::ConstElemIterator end(env.elements_end());
      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;
        if (_quantum_optics->includes_region(elem->subdomain_id()))
        {

        }
      }
    }
  }
}
