// $Id: MaxwellAbsorption.C 2063 2010-09-03 13:11:49Z maufder $

#include "MaxwellAbsorption.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include <assert.h>

#include "equation_systems.h"
#include "limits.h"
#include "Utils.h"
#include "point.h"
#include "math.h"
#include "Constants.h"
#include "SimulationInterface.h"

#include "TiberModule.h"


MaxwellAbsorption::MaxwellAbsorption(const ModelOptions& options)
 : SimulationInterface(options)
{
  has_solution_vector(false);
  increment_solve_sequence_number(); //This simulation is solved by default making get_solution() accessable.
}


void MaxwellAbsorption::do_init() {
  std::string m = get_options().get_option("maxwell_simulation", "maxwell_boundary");
  maxwell = SimulationInterface::find_simulation(m);

  std::string r_sim = get_options().get_option("referenced_points_simulation", "excitontransport");
  referencedPointsSimulation = SimulationInterface::find_simulation(r_sim);

  lambdaStart = get_options().get_option("lambdaStart", 400);
  lambdaEnd = get_options().get_option("lambdaEnd", 700);
  lambdaStep = get_options().get_option("lambdaEnd", 10);
}


//=======================================================================================================//
void MaxwellAbsorption::do_solve() {
  const MeshBase& mesh = get_mesh();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_local_elements_end();

  ID id = maxwell->get_solution_id("Efield");

  for (double lambda = lambdaStart; lambda < lambdaEnd; lambda++) {
    double W = 2 * M_PI / lambda * Constants::c;
    double W_in_eV = Constants::hbar * W / Constants::e;

    maxwell->get_options().set_option("W", W_in_eV);
    // set power (?) how

    maxwell->solve();

    MeshBase::const_element_iterator el =
                                    mesh.active_local_elements_begin();

    // loop over all active elements
    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;

      std::map<ID, std::vector<double>> maxwell_solutions;
      // add sol vars to maxwell_solutions

      maxwell->get_solution(elem, maxwell_solutions, integration_points[elem->centroid()]);

      for (int i = 0; i < maxwell_solutions[id].size(); i++) {
          //solution[integration_points[elem->centroid()][i]] += ; //smth like that
      }
    }
  }
}


void
MaxwellAbsorption::do_setup_solution_variables(void)
{
  declare_solution_ext("Absorption", Absorption, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
}

void MaxwellAbsorption::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<Point>& points)
{
  if (solutions.count(Absorption)) {
    std::vector<double>& sol = solutions[Absorption];
    sol.resize(points.size());

    for (int i = 0; i < points.size(); i++) {
      if (!solution.count(points[i])) {
        integration_points[elem->centroid()].push_back(points[i]);
        solution[points[i]] = 0;
      }

      sol[i] = solution[points[i]];
    }
  }
}
