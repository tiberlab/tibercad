// $Id$

#include "PoissonDirichlet.h"

#include "TiberModule.h"
#include "SimulationInterface.h"

using namespace std;
using namespace libMesh;


void
PoissonDirichlet::do_init(void)
{
  string str = get_option("potential", "0");
  istringstream is(str);

  double val;
  if ((is >> val) || (str[0] == '$'))
  {
     get_parameter("potential", _potential);
  }
  else
  {
     vector<string> pot;
     Utils::extract_vector(str, pot);
     _sim.resize(1);
     _sim[0] = SimulationInterface::find_solution_provider(pot[0],"ElPotential");

     if (_sim[0].first == NULL || _sim[0].second == INVALID_ID)
     {
         throw InitFailedException(pot[0] + " is invalid identifier for "
                "a module providing ElPotential");
     }
  }
}


void
PoissonDirichlet::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  if (_sim.size() > 0)
  {
     ID id = _sim[0].second;
     SimulationInterface* si = _sim[0].first;
     si->get_solution(elem, id, _potential, point);
  }
  
  set_coefficients(1, 0, _potential);
}
