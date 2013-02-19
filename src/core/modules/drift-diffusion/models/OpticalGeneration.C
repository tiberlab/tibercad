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

  string gen_str(get_option("generation", "0"));
  istringstream is(gen_str);
//  string gen_file(get_option("gen_file", "0"));
//  istringstream is2(gen_file);
 
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
    recomb_e = recomb_h = -_generation;
}




double
//OpticalGeneration::read_file(char *filename, double p)
OpticalGeneration::read_file(void)
{
 //This simplified generation model reads the generation from an external file.
 //The structure of the file is (x, y, z, gen) and the coordinates of the nodes must be
 //the same as for the element.
/*
  
  bool coord = false;
  double x1 = 0.0;
  double x2 = 0.0;
  double x3 = 0.0;
// Qui mettere le coordinate vere dal modello

  string str;
  char cstr[200];
  char* tok;
  int column, i, j;
*/
  double gen;
/*
  vector<double> x;

  ifstream file;
  file.open (filename);

  //scroll file up to the end # lines
  while (! file.eof() )
  {
    getline(file, str);
    strcpy(cstr, str.c_str());
    if((int)cstr[0] != 35)
    {
      break;
    }
  }

  //this counts the number of tokens
  tok = strtok (cstr, " ");
  i = 0;
  while (tok != NULL)
  {
    i++;
    tok = strtok (NULL, " ");
  }
  //define number of columns
  column = i;
  x.resize(column);
  
  //first line is already buffered so just insert in vectors, then continue to read the stream
  i=0;
  strcpy(cstr, str.c_str());
  tok = strtok(cstr," ");

  while (tok != NULL)
  {
    x[i].push_back(atof(tok));
    tok = strtok (NULL, " ");
    i ++;
  }

  if( (x[0] == x1) && (x[1] == x2) && (x[2] == x3) )
  {
    coord = true;
    gen = x[3];
  } 

  while( (!file.eof()) && (coord == false))
  {
    getline(file, str);
    strcpy(cstr, str.c_str());
    if (!isprint(cstr[0])) {
      break;
    }
    tok = strtok(cstr," ");
    i = 0;
    while (tok != NULL)
    {
      x[i].push_back(atof(tok));
      tok = strtok(NULL," ");
      i++;
    }
    
    if( (x[0] == x1) && (x[1] == x2) && (x[2] == x3) )
    {
      coord = true;
      gen = x[3];
    } 
  }
*/
  return gen;
}

