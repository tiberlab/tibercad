// $Id$


#include "OpticalGeneration.h"
#include "Material.h"
#include "SimulationInterface.h"
#include "DriftDiffusionProperties.h"

#include <string>
#include "elem.h"


TIBER_MODULE(OpticalGeneration, recombination, optical)

using namespace std;



void
OpticalGeneration::do_init(void)
{
  // G is a sweepable value, so check it!
  //get_parameter("G", G_);

  string gen_str(get_option("generation", ""));
  istringstream is(gen_str);
  double val;
  if ((is >> val) || (gen_str[0] == '$'))
  {
    get_parameter("generation", _generation);
  }
  else
  {
    vector<string> gens;
    Utils::extract_vector(gen_str, gens);
    _generation_model.resize(gens.size());
    _gen_id.resize(gens.size());
    for (size_t i = 0; i < gens.size(); ++i)
    {
      _generation_model[i] = SimulationInterface::find_simulation(gens[i]);
      if (_generation_model[i] == NULL)
        throw InitFailedException("Cannot find generation model: " + gens[i]);

      _gen_id[i] = _generation_model[i]->get_solution_id("Generation");

    }
   }
}



void
OpticalGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  recomb_e[0] = recomb_h[0] = 0;
  recomb_e[1] = recomb_h[1] = 0;
}

void
OpticalGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
    if (_generation_model.size() > 0)
    {
      
      DriftDiffusionProperties& dd = get_driftdiffusionproperties();

      const Elem* el = dd.get_element();

      _generation = 0;

      vector<double> tmp(1);
      for (size_t i = 0; i < _generation_model.size(); ++i)
      {
        if (_generation_model[i]->get_solution(el, _gen_id[i], tmp,
            vector<Point>(1, dd.get_coordinates())))
          _generation += tmp[0];
      }

    }

    recomb_e = recomb_h = -_generation;
}

