#include "ExcitonTransport.h"
#include "ExcitonDissociation.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"



ExcitonDissociation::ExcitonDissociation(void)
  : RecombinationModelInterface(),
    _d(1.0),
    _a(0.0),
    _exciton_sim(NULL)

{
}

void
ExcitonDissociation::set_model_options(const ModelOptions& options)
{
  ModelOptions::const_iterator it = options.find("damping");
  if (it != options.end())
    _d = atof((it->second).c_str());
  
  it = options.find("trapping_probability");
  if (it != options.end())
    _a = atof((it->second).c_str());
}

void
ExcitonDissociation::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  recomb_e = 0.0;
  
  if (_exciton_sim != NULL)
  {
    const Elem* el = dd.get_element();

    assert(_exciton_sim->get_exciton_model());
    _exciton_sim->get_exciton_model()->reinit(el, &dd);

    const Point& p = *(dd.get_coordinates());
    double u = _exciton_sim->get_solution(el, p);
    _exciton_sim->get_exciton_model()->calculate_densities(u);
    //_exciton_sim->get_exciton_model()->calculate_recombination_rate();
    //recomb_e =
    //  -_d * _exciton_sim->get_exciton_model()->get_recombination_rate();
    recomb_e = -_d * (1 - _a) *
      _exciton_sim->get_exciton_model()->get_nonradiative_recombination_rate();
  }
  
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

const std::string
ExcitonDissociation::get_name(void) const
{
  return "exciton_dissociation";
}
