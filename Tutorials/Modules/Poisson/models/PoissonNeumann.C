// $Id: PoissonNeumann.C 2029 2010-07-21 13:01:58Z maufder $

#include "PoissonNeumann.h"

#include "TiberModule.h"
#include "SimulationInterface.h"
#include "mesh.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"

using namespace std;
using namespace libMesh;

void
PoissonNeumann::do_init(void)
{
  
  string str = get_option("displacement", "0");
  istringstream is(str);

  double val;
  if ((is >> val) || (str[0] == '$'))
  {
     // taken in C/m^2
     get_parameter("displacement", _field);
  }
  else
  {
     vector<string> displ;
     Utils::extract_vector(str, displ);
     _sim.resize(1);
     _sim[0] = SimulationInterface::find_solution_provider(displ[0],"Displacement");

     if (_sim[0].first == NULL || _sim[0].second == INVALID_ID)
     {
         throw InitFailedException(displ[0] + " is invalid identifier for "
                "a module providing Displacement vector");
     }

     dim = (SimulationInterface::get_simulation(get_simulator_id())->get_mesh()).mesh_dimension(); 

  }
}


void
PoissonNeumann::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  if (_sim.size() > 0)
  {
     vector<double> displ(3,0.0);
     ID id = _sim[0].second;
     SimulationInterface* si = _sim[0].first;
     si->get_solution(elem, id, displ, point);
     //construct normal component
     unique_ptr<FEBase> fe(FEBase::build(dim, FEType()));
     QGauss qrule(dim-1, CONSTANT); // Order 0 rule because in this way we take centroid's normal
     fe->attach_quadrature_rule(&qrule);
     const std::vector<Point>& normal = fe->get_normals();
     fe->reinit(elem, side);
     _field = -normal[0](0)*displ[0]-normal[0](1)*displ[1]-normal[0](2)*displ[2];     
  }
  
  set_coefficients(0, 1, _field);

}
