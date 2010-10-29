// $Id$

#include "StrainInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"
#include "RuntimeException.h"
#include "Messages.h"

#include "elem.h"

#include "tensor.h"



StrainInterface::StrainInterface(void)
  : _simulation(NULL)
{
}



bool
StrainInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _strain_id = _simulation->get_solution_id("Strain");
    _strain_cryst_id = _simulation->get_solution_id("StrainCrystal");
    _stress_id = _simulation->get_solution_id("Stress");
    _stress_cryst_id = _simulation->get_solution_id("StressCrystal");

    if (_strain_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'Strain\'.");

    if (_strain_cryst_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'StrainCrystal\'.");

    if (_stress_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'Stress\'.");

    if (_stress_cryst_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'StressCrystal\'.");

    
    for (int i = 0; i < 9; i++)
      if ((_strain_id == INVALID_ID) && (_strain_cryst_id == INVALID_ID) &&
          (_stress_id == INVALID_ID) && (_stress_cryst_id == INVALID_ID))
        throw InitFailedException("Simulation " + name +
            " is missing strain related variables");

    answer = true;
  }

  return answer;
}



void
StrainInterface::get_strain(const Elem* elem, const Point& point, Tensor2Sym& strain)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'Strain\'.");

    _get_data(elem, point, strain, _strain_id);
  }
}


void
StrainInterface::get_crystal_strain(const Elem* elem, const Point& point, Tensor2Sym& strain)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'StrainCrystal\'.");

    _get_data(elem, point, strain, _strain_cryst_id);
  }
}


void
StrainInterface::get_stress(const Elem* elem, const Point& point, Tensor2Sym& stress)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'Stress\'.");

    _get_data(elem, point, stress, _stress_id);
  }
}


void
StrainInterface::get_crystal_stress(const Elem* elem, const Point& point, Tensor2Sym& stress)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'StressCrystal\'.");

    _get_data(elem, point, stress, _stress_cryst_id);
  }
}



void
StrainInterface::_get_data(const Elem* elem, const Point& point,
    Tensor2Sym& data, ID id)
{
  std::vector<Point> p(1);
  p[0] = point;
  std::vector<double> values(6);

  if (_simulation->get_solution(elem, id, values, p))
  {
    data(1,1) = values[0];
    data(2,2) = values[1];
    data(3,3) = values[2];
    data(2,1) = values[3];
    data(3,2) = values[4];
    data(3,1) = values[5];
  }
}


