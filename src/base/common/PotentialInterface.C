// $Id$

#include "PotentialInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include "elem.h"



PotentialInterface::PotentialInterface(void) :
  _simulation(NULL),
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

    _id = _simulation->get_solution_id("ElPotential");

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable 'ElPotential'");


    //TODO: THIS PART IS ADDED FOR DIRTY IWCE CALCULATIONS!!!
    //------------------------------------------------------------
    _id_chem_el = _simulation->get_solution_id("eQFermi");

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable 'eQFermi'");


    _id_chem_hl = _simulation->get_solution_id("hQFermi");

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable 'hQFermi'");

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

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id, potentials))
  {
    int nn = elem->n_nodes();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


void
PotentialInterface::get_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials,
    bool local_coord)
{
  assert(elem != NULL);

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id, potentials, p, local_coord))
  {
    int nn = p.size();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


double
PotentialInterface::get_potential(const Elem* elem, const Point& p,
    bool local_coord)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_potential(elem, ps, temp, local_coord);

  return temp[0];
}


void
PotentialInterface::get_el_chem_potential(const Elem* elem,
    std::vector<double>& potentials)
{
  assert(elem != NULL);

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id_chem_el, potentials))
  {
    int nn = elem->n_nodes();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


void
PotentialInterface::get_el_chem_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials,
    bool local_coord)
{
  assert(elem != NULL);

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id_chem_el, potentials, p, local_coord))
  {
    int nn = p.size();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


double
PotentialInterface::get_el_chem_potential(const Elem* elem, const Point& p,
    bool local_coord)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_el_chem_potential(elem, ps, temp, local_coord);

  return temp[0];
}


void
PotentialInterface::get_hl_chem_potential(const Elem* elem,
    std::vector<double>& potentials)
{
  assert(elem != NULL);

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id_chem_hl, potentials))
  {
    int nn = elem->n_nodes();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


void
PotentialInterface::get_hl_chem_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potentials,
    bool local_coord)
{
  assert(elem != NULL);

  if ((_simulation == NULL) ||
      !_simulation->get_solution(elem, _id_chem_hl, potentials, p, local_coord))
  {
    int nn = p.size();
    potentials.resize(nn);
    for (int i = 0; i < nn; i++)
      potentials[i] = 0.0;
  }
}


double
PotentialInterface::get_hl_chem_potential(const Elem* elem, const Point& p,
    bool local_coord)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_hl_chem_potential(elem, ps, temp, local_coord);

  return temp[0];
}
