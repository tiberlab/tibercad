// $Id$

#include "DSSCGeneration.h"
#include "TiberLinearSystem.h"

#include "mesh_base.h"


TIBER_MODULE(DSSCGeneration, dssc_generation)


DSSCGeneration::DSSCGeneration(const ModelOptions& options) :
  SimulationInterface(options),
  _direction(0),
  _intensity(0)
{
  int dim = get_mesh().mesh_dimension();
  switch (dim)
  {
    case 2:
      _direction(1) = -1;
      break;

    case 3:
      _direction(2) = -1;
      break;

    default:
      _direction(0) = 1;
      break;
  }

}



void
DSSCGeneration::parse_options(void)
{
}


void
DSSCGeneration::do_init(void)
{

  get_parameter("light_direction", _direction);
  get_parameter("light_intensity", _intensity);


  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  //system.attach_assembly_routine(assemble_system);


  system.add_variable("d", libMeshEnums::FIRST);

  system.add_vector("G");

  system.init();
}



void
DSSCGeneration::do_setup_solution_variables(void)
{
  declare_solution(Generation, REAL, NODES, "1/(cm^3*s)");
}



void
DSSCGeneration::do_solve(void)
{

  if (!_d_calculated)
    _calculate_distances();



}


void
DSSCGeneration::_calculate_distances(void)
{

  _d_calculated = true;
}

