// $Id$

#include "TemperatureInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include "elem.h"

std::string
TemperatureInterface::_variable_name = "temperature";



TemperatureInterface::TemperatureInterface(void)
  : _simulation(NULL),
    _id(INVALID_ID)
{
}



void
TemperatureInterface::set_simulation(const std::string& name)
{
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _id_set.clear();
    _id = _simulation->get_variable_id(_variable_name);

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable '" + _variable_name + "'");

    _id_set.insert(_id);
  }
}


void
TemperatureInterface::get_temperature(const Elem* elem,
    std::vector<double>& temperatures)
{
  assert(elem != NULL);
  
  int nn = elem->n_nodes();

  temperatures.resize(nn);

  if (_simulation == NULL)
  {
    for (int i = 0; i < nn; i++)
      temperatures[i] = SimulationOptions::temperature;
  }
  else
  {
    
    std::vector<std::map<ID, double> > temp;
    if (_simulation->get_solution(elem, _id_set, temp))
      for (int i = 0; i < nn; i++)
        temperatures[i] = temp[i][_id];
    else
      for (int i = 0; i < nn; i++)
        temperatures[i] = SimulationOptions::temperature;
  }
}


void
TemperatureInterface::get_temperature(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& temperatures)
{
  assert(elem != NULL);

  int nn = p.size();

  temperatures.resize(nn);

  if (_simulation == NULL)
  {
    for (int i = 0; i < nn; i++)
      temperatures[i] = SimulationOptions::temperature;
  }
  else
  {
    std::vector<std::map<ID, double> > temp;
    if (_simulation->get_solution(elem, p, _id_set, temp))
      for (int i = 0; i < nn; i++)
        temperatures[i] = temp[i][_id];
    else
      for (int i = 0; i < nn; i++)
        temperatures[i] = SimulationOptions::temperature;
  }
}


double
TemperatureInterface::get_temperature(const Elem* elem, const Point& p)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_temperature(elem, ps, temp);
  
  return temp[0];
}
