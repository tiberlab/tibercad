
#include "KmcInterface.h"
#include "SimulationEnvironment.h"
#include "SolutionDescriptor.h"
#include "Messages.h"


void
KmcInterface::do_init()
{

}



bool
KmcInterface::get_boundary_potentials(Boundary* bd, 
                                      const Point& pt,
                                      double& phi,
                                      double& efermi,
                                      double& hfermi)
{
  phi = 0.0;
  efermi = 0.0;
  hfermi = 0.0;

  if (_pot_sol.first == NULL || _mue_sol.first == NULL || _muh_sol.first == NULL) return false;

  SimulationInterface *pot_sim = _pot_sol.first;
  SimulationInterface *mue_sim = _mue_sol.first;
  SimulationInterface *muh_sim = _muh_sol.first;

  // Problem 1: the side I get from SimulationEnvironment are on the 
  SimulationEnvironment& env = get_environment();
  SimulationEnvironment::BoundarySideIterator it = env.boundary_sides_begin(bd->get_name());
  SimulationEnvironment::BoundarySideIterator itend = env.boundary_sides_end(bd->get_name());

  const Elem* neighbour = ((*it).elem())->neighbor((*it).side());

  if (neighbour == NULL)
  {
      Messages::warning((bd->get_name())+" is not an inner boundary");
      return false;
  }

  // Iterate on all element side of a given boundary :o !
  for (; it != itend; ++it)
  {

     const ElementSide& side = *it;
     const Elem* neigh = side.elem()->neighbor(side.side());  
     double val;

     if (side.elem()->contains_point(pt) || neigh->contains_point(pt))
     {
        pot_sim->get_solution(neigh, _pot_sol.second, val, neigh->centroid());
        phi = val;
        mue_sim->get_solution(neigh, _mue_sol.second, val, neigh->centroid());
        efermi = val;
        muh_sim->get_solution(neigh, _muh_sol.second, val, neigh->centroid());
        hfermi = val;
     } 
 
  }

  return true;

}


