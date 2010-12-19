// $Id$

#include "TemperatureInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include "elem.h"

std::string
TemperatureInterface::_variable_name = "LatticeTemp";



TemperatureInterface::TemperatureInterface(void)
  : _simulation(NULL),
    _id(INVALID_ID)
{
}



bool
TemperatureInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _id = _simulation->get_solution_id(_variable_name);

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable '" + _variable_name + "'");

    answer = true;
  }

  return answer;
}


void
TemperatureInterface::get_temperature(const Elem* elem,
    std::vector<double>& temperatures)
{
  assert(elem != NULL);

  int nn = elem->n_nodes();
  std::vector<Point> nodes(nn);

  for (unsigned int i = 0; i < nn; ++i)
    nodes[i] = elem->point(i);

  get_temperature(elem, nodes, temperatures);
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
    std::map<ID, std::vector<double> > temp;
    temp[_id] = std::vector<double>();

    if (_simulation->get_solution(elem, temp, p))
      temperatures = temp[_id];
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

bool
TemperatureInterface::is_solved(void) const
{
  return (_simulation->is_solved());
}
