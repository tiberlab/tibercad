
#include "OpticalRecombination.h"
#include "TiberModule.h"
#include "SimulationInterface.h"
#include "Messages.h"
#include <string>



void 
OpticalRecombination::calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda)
{
  string recombination_str(get_option("radiative_recombination", "0"));
   vector<string> recombination;
   Utils::extract_vector(recombination_str, recombination);
   _recombination_model.resize(recombination.size());
   _recombination_id.resize(recombination.size());
   for (size_t i = 0; i < recombination.size(); ++i)
    {
      pair<SimulationInterface*, ID> provider =
        SimulationInterface::find_solution_provider(recombination[i]);
      _recombination_model[i] = provider.first;
      _recombination_id[i] = provider.second;
      if (_recombination_model[i] == NULL)
        throw InitFailedException("Cannot find radiative_recombination model: " + recombination[i]);
    }
	
    _recombination_rate = 0;
	

    vector<double> tmp(1);
    for (size_t i = 0; i < _recombination_model.size(); ++i)
    {
      if (_recombination_model[i]->get_solution(elem, _recombination_id[i], tmp,
          vector<Point>(1, point)))
        _recombination_rate += tmp[0];
    }

	

	set_emission_power(_recombination_rate ); 
}

