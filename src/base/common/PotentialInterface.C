// $Id$

#include "PotentialInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "Atom.h"
#include "InitFailedException.h"

#include "elem.h"
#include "mesh_base.h"


PotentialInterface::PotentialInterface(void) :
  _simulation(NULL),
  _id(INVALID_ID)
{
}


PotentialInterface::PotentialInterface(const std::string& name, const std::string& variable) :
  _simulation(NULL),
  _id(INVALID_ID)
{
  set_simulation(name, variable); 
}



bool
PotentialInterface::set_simulation(const std::string& name, const std::string& variable)
{
  bool answer = false;
  if (name != "")
  {
    SimulationInterface::SolutionProvider prov =
        SimulationInterface::find_solution_provider(name, variable);

    _simulation = prov.first;
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _id = prov.second;

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable '"+variable+"'");

    answer = true;
  }

  return answer;
}



void
PotentialInterface::get_potential(const libMesh::Elem* elem,
    const std::vector<libMesh::Point>& p, std::vector<double>& potentials,
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
PotentialInterface::get_potential(const libMesh::Elem* elem, const libMesh::Point& p,
    bool local_coord)
{
  std::vector<libMesh::Point> ps(1, p);
  std::vector<double> temp(1);

  get_potential(elem, ps, temp, local_coord);

  return temp[0];
}



double
PotentialInterface::get_potential(const Atom* atom)
{
  libMesh::Point p(0);
  switch (get_simulation()->get_mesh().mesh_dimension())
  {
    case 3:
      p(2) = atom->get_position(2);
    case 2:
      p(1) = atom->get_position(1);
    default:
      p(0) = atom->get_position(0);
  }
  return get_potential(atom->get_elem(), p);
}


