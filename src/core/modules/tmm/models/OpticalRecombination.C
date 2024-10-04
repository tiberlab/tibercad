
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
	get_parameter("emission_peak", _emission_wavelength);
	get_parameter("multiplier", _multiplier);
	
	if (_multiplier > 1 || _multiplier<= 0)
		_multiplier = 1;
	
	
	//std::cout << "_multiplier " << _multiplier << endl;
	  
	double plank_const = 6.62607015e-34;  //Planck constant[J*s]
	double c0 = 2.998e8 * 1e9;  //speed of ligth [nm/s]
	double power;
	//energy of each photon multiplyed by number of recombinations
	power = _recombination_rate ;
	set_emission_power(power ); 
}

