#include "GuestInterface.h"
#include "DriftDiffusionProperties.h"
#include "TiberModule.h"
#include "SimulationInterface.h"


GuestInterface::GuestInterface(const ModelOptions& options)
  : ElectricalContact(options)
{
  has_current(false);
  set_type(0, NEUMANN);
  set_type(1, NEUMANN);
  set_type(2, NEUMANN);
}


void
GuestInterface::do_init(void)
{
  ElectricalContact::do_init();

  std::string sim_name = get_option("host_simulation", "");


  if (sim_name != "") {
    SimulationInterface* sim = SimulationInterface::find_simulation(sim_name);

    if (sim == NULL) {
      std::string msg("External host simulation " + std::string(sim_name) + " not found");
      throw InitFailedException(msg);
    }
    _host_sim = sim->get_id();
    set_type(0, DIRICHLET);
  }
}


void
GuestInterface::do_compute(void)
{

  if (_host_sim > 0)
  {
    SimulationInterface* sim = SimulationInterface::get_simulation(_host_sim);
    if (sim->is_solved())
    {
      DriftDiffusionProperties* mod = sim->get_bulk_model<DriftDiffusionProperties>(get_element());
      if ( mod != NULL )
      {
        double Ech  = mod->get_conduction_band().get_band_edge();
        ID Ec_id = sim->get_solution_id("Ec");

        double x = 0.0;
        sim->get_solution(get_element(), Ec_id, x, get_coordinates());
        x = Ech - x;
        coeff_g(0) = x;
      }
    }
  }

}