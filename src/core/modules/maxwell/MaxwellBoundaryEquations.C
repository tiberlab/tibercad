// $Id: MaxwellBoundaryEquations.C 2063 2010-09-03 13:11:49Z maufder $

#include "MaxwellBoundaryEquations.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "MaxwellBoundaryProperties.h"

#include "equation_systems.h"
#include "dense_submatrix.h"
#include "limits.h"
#include "Utils.h"
#include "TiberLinearSystem.h"
#include "VectorFEBase1D.h"

using namespace std;
using namespace Constants;

#include "sys/time.h"

TIBER_MODULE(MaxwellBoundaryEquations, maxwell, boundary)

MaxwellBoundaryEquations::MaxwellBoundaryEquations(const ModelOptions& options)
 : SimulationInterface(options)
{
  has_solution_vector(false);
}

//=====================================================//
BoundaryProperties* MaxwellBoundaryEquations::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  return MaxwellBoundaryProperties::create(options);
}

//=======================================================================================================//
PhysicalModel*  MaxwellBoundaryEquations::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  MaxwellPhysicalModel* model = dynamic_cast<MaxwellPhysicalModel*> ( PhysicalModelInterface::create("maxwell", mat, options) );

  if (model == NULL)
    throw ModelErrorException("MaxwellBoundaryEquations: cannot create MaxwellPhysicalModel");

  return(model);

}

//=======================================================================================================//

void MaxwellBoundaryEquations::do_init() {
  W = get_options().get_option("W", 0.0) / Constants::hbar * Constants::e;

  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  system.add_variable("Ex", FIRST);
  system.add_variable("Ey", FIRST);
  system.add_variable("Ez", FIRST);

  system.attach_assemble_function(assemble_maxwell_equations);
  system.init();
}


//=======================================================================================================//
void MaxwellBoundaryEquations::do_solve() {
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  //system.solve();
  VectorFEBase1D test(1.0);
  std::cout << "TEST FOUR" << test.FOUR << "\n";
}

void
MaxwellBoundaryEquations::do_setup_solution_variables(void)
{
  declare_solution(Efield, REAL, NODES, "abs");
}

void
MaxwellBoundaryEquations::assemble_maxwell_equations(EquationSystems& es,
    const std::string& system_name)
{


/*  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();*/

}

void MaxwellBoundaryEquations::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<Point>& points)
{

}

