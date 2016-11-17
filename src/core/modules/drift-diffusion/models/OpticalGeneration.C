// $Id$


#include "OpticalGeneration.h"
#include "Material.h"
#include "SimulationInterface.h"
#include "DriftDiffusionProperties.h"
#include "Messages.h"

#include "TiberModule.h"

#include <string>
#include "elem.h"



using namespace std;

void
OpticalGeneration::do_init(void)
{
  // G is a sweepable value, so check it!
  //get_parameter("G", G_);

  string gen_str(get_option("generation", "0"));
  istringstream is(gen_str);

  get_parameter("multiplier", _multiplier);
 
  double val;
  if ((is >> val) || (gen_str[0] == '$'))
  {
  // here we define if the generation is directly fixed by a number or 
  // must be read in a data file. If gen_file != 0 then the program looks for a file.
  // In this latter case generation contains the number of sun. The total generation
  // is equal to _generation = _sun * G (from the file).

//    if (gen_file[0] != '0')
//    {
//      _read_file = true;
//      get_parameter("generation", _sun); 
//    }
//    else
//    {
      get_parameter("generation", _generation);
//    }
  }
  else
  {
    vector<string> gens;
    Utils::extract_vector(gen_str, gens);
    _generation_model.resize(gens.size());
    _gen_id.resize(gens.size());
    for (size_t i = 0; i < gens.size(); ++i)
    {
      pair<SimulationInterface*, ID> provider =
        SimulationInterface::find_solution_provider(gens[i]);
      _generation_model[i] = provider.first;
      _gen_id[i] = provider.second;
      cerr <<_generation_model[i] << " " << _gen_id[i] << endl;

      if (_generation_model[i] == NULL)
        throw InitFailedException("Cannot find generation model: " + gens[i]);

      if (_gen_id[i] == INVALID_ID)
        _gen_id[i] = _generation_model[i]->get_solution_id("generation");
      if (_gen_id[i] == INVALID_ID)
        _gen_id[i] = _generation_model[i]->get_solution_id("Generation");
      if (_gen_id[i] == INVALID_ID)
        throw InitFailedException("Cannot find a solution variable for generation model: " + gens[i]);


      ostringstream os;
      os << "Found generation model: " << gens[i] << " (" << _generation_model[i] << "), variable ID " << _gen_id[i] << endl;
      Messages::info(os.str());

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
/*
    if (_read_file == true)
    {
      DriftDiffusionProperties& dd = get_driftdiffusionproperties();
      const Elem* el = dd.get_element();
      vector<Point> p = (1,dd.get_coordinates());
      // prendere le cooridnate e passarle al file read_file
      _generation = _sun * read_file(&gen_file, 1);
    }
*/
  }

  recomb_e = recomb_h = -_multiplier * _generation;
}



