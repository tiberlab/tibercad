// $Id$

#include "PotentialInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include "elem.h"

std::string
PotentialInterface::_variable_name = "ElPotential";



PotentialInterface::PotentialInterface(void)
: _simulation(NULL),
_id(INVALID_ID)
{
}



bool
PotentialInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  if (name != "")
    {
      _simulation = SimulationInterface::find_simulation(name);
      if (_simulation == NULL)
        throw InitFailedException("No such simulation found: " + name);

      _id_set.clear();
      _id = _simulation->get_solution_id(_variable_name);

      if (_id == INVALID_ID)
        throw InitFailedException("Simulation " + name +
            " has no variable '" + _variable_name + "'");

      _id_set.insert(_id);

      //TODO: THIS PART IS ADDED FOR DIRTY IWCE CALCULATIONS!!!
      //------------------------------------------------------------
      _id_chem_el = _simulation->get_solution_id("QFermi_e");

      if (_id == INVALID_ID)
        throw InitFailedException("Simulation " + name +
            " has no variable '" + _variable_name + "'");

      _id_set.insert(_id_chem_el);

      _id_chem_hl = _simulation->get_solution_id("QFermi_h");

      if (_id == INVALID_ID)
        throw InitFailedException("Simulation " + name +
            " has no variable '" + _variable_name + "'");

      _id_set.insert(_id_chem_hl);
      //------------------------------------------------------------


      answer = true;
    }

  return answer;
}


void
PotentialInterface::get_potential(const Elem* elem,
    std::vector<double>& potentials)
{
  assert(elem != NULL);

  int nn = elem->n_nodes();

  potentials.resize(nn);


  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {

      std::vector<std::map<ID, double> > temp;

      if (_simulation->get_solution(elem, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id];

      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


void
PotentialInterface::get_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials)
{
  assert(elem != NULL);

  int nn = p.size();

  potentials.resize(nn);

  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {
      std::vector<std::map<ID, double> > temp;
      if (_simulation->get_solution(elem, p, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id];
      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


double
PotentialInterface::get_potential(const Elem* elem, const Point& p)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_potential(elem, ps, temp);

  return temp[0];
}


void
PotentialInterface::get_el_chem_potential(const Elem* elem,
    std::vector<double>& potentials)
{
  assert(elem != NULL);

  int nn = elem->n_nodes();

  potentials.resize(nn);


  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {

      std::vector<std::map<ID, double> > temp;

      if (_simulation->get_solution(elem, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id_chem_el];

      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


void
PotentialInterface::get_el_chem_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials)
{

  assert(elem != NULL);

  int nn = p.size();

  potentials.resize(nn);

  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {
      std::vector<std::map<ID, double> > temp;
      if (_simulation->get_solution(elem, p, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id_chem_el];
      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


double
PotentialInterface::get_el_chem_potential(const Elem* elem, const Point& p)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_el_chem_potential(elem, ps, temp);

  return temp[0];
}


void
PotentialInterface::get_hl_chem_potential(const Elem* elem,
    std::vector<double>& potentials)
{
  assert(elem != NULL);

  int nn = elem->n_nodes();

  potentials.resize(nn);


  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {

      std::vector<std::map<ID, double> > temp;

      if (_simulation->get_solution(elem, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id_chem_hl];

      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


void
PotentialInterface::get_hl_chem_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials)
{

  assert(elem != NULL);

  int nn = p.size();

  potentials.resize(nn);

  if (_simulation == NULL)
    {
      for (int i = 0; i < nn; i++)
        potentials[i] = 0.0;
    }
  else
    {
      std::vector<std::map<ID, double> > temp;
      if (_simulation->get_solution(elem, p, _id_set, temp))
        for (int i = 0; i < nn; i++)
          potentials[i] = temp[i][_id_chem_hl];
      else
        for (int i = 0; i < nn; i++)
          potentials[i] = 0.0;
    }
}


double
PotentialInterface::get_hl_chem_potential(const Elem* elem, const Point& p)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_hl_chem_potential(elem, ps, temp);

  return temp[0];
}
