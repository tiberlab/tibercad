// $Id$

#include "ExcitonDissociation.h"

#include "ExcitonTransport.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Utils.h"


TIBER_MODULE(ExcitonDissociation, exciton_dissociation)



void
ExcitonDissociation::do_init(void)
{
  d_ = get_options().get_option("damping", 1.0);

  std::string ex = get_options().get_option("exciton_simulation",
      Utils::extract_typename(typeid(exciton_sim_)));

  // find the exciton simulation to use
  exciton_sim_ = dynamic_cast<ExcitonTransport*>(
      SimulationInterface::find_simulation(ex));

  if (exciton_sim_ == NULL)
  {
    std::string msg("ExcitonDissociation: Simulation " +
        std::string(ex) + " not found");
    throw InitFailedException(msg);
  }

  _Rdiss_id = exciton_sim_->get_variable_id("ChemPot");
}

void
ExcitonDissociation::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  const Elem* el = dd.get_element();

  // we only use the exciton simulation if it has been solved before
  if (exciton_sim_->is_solved())
  {
    ID ex_id = exciton_sim_->get_id();
    ExcitonProperties* mod =
      static_cast<ExcitonProperties*>(get_material()->get_model(ex_id));
    mod->reinit(el);

    double x = 0.0;
    bool succ = exciton_sim_->get_solution(el, dd.get_coordinates(), _Rdiss_id, x);
    if (succ)
      recomb_e = -d_ * mod->get_dissociation_rate();
  }
  else
    recomb_e = 0.0;

  recomb_h = recomb_e;
}

void
ExcitonDissociation::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
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

