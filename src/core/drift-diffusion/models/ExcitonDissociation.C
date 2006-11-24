// $Id$

#include "ExcitonTransport.h"
#include "ExcitonDissociation.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Utils.h"


void
ExcitonDissociation::do_init(void)
{
  _d = get_options().get_option("damping", 1.0);

  std::string ex = get_options().get_option("exciton_simulation",
      Utils::extract_typename(typeid(_exciton_sim)));

  // find the exciton simulation to use
  _exciton_sim = dynamic_cast<ExcitonTransport*>(
      SimulationInterface::find_simulation(ex));

  if (_exciton_sim == NULL)
  {
    std::string msg("Simulation "+std::string(ex)+" not found");
    throw InitFailedException(msg);
  }
}

void
ExcitonDissociation::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  const Elem* el = dd.get_element();

  // we only use the exciton simulation if it has been solved before
  if (_exciton_sim->is_solved())
  {
    ID ex_id = _exciton_sim->get_id();
    ExcitonProperties* mod =
      static_cast<ExcitonProperties*>(get_material()->get_model(ex_id));
    mod->reinit(el);

    double u = _exciton_sim->get_solution(el, dd.get_coordinates());
    mod->set_effective_potential(u);
    mod->calculate_density();
    recomb_e = -_d * mod->get_dissociation_rate();
  }
  else
    recomb_e = 0.0;

  recomb_h = recomb_e;
}

void
ExcitonDissociation::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  recomb_e[0] = recomb_h[0] = 0.0;
  recomb_e[1] = recomb_h[1] = 0.0;
  recomb_e[2] = recomb_h[2] = 0.0;
}


void
ExcitonDissociation::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  ignore_unused_variable(comp_A);
  ignore_unused_variable(comp_B);
  ignore_unused_variable(xa);
}

