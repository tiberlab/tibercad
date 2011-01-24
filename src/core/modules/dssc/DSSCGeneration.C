// $Id$

#include "DSSCGeneration.h"
#include "TiberLinearSystem.h"
#include "SimulationEnvironment.h"

#include "mesh_base.h"
#include "dof_map.h"


TIBER_MODULE(DSSCGeneration, dssc_generation)


using namespace std;


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

DSSCGeneration::~DSSCGeneration(void)
{

}


void
DSSCGeneration::parse_options(void)
{
}


void
DSSCGeneration::do_init(void)
{

  get_parameter("light_direction", _direction);
  double len = _direction.size();
  _direction /= len;

  get_parameter("light_intensity", _intensity);


  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  //system.attach_assembly_routine(assemble_system);


  system.add_variable("d", libMeshEnums::FIRST);

  system.add_vector("G");

  system.init();

}

PhysicalModel*
DSSCGeneration::create_bulk_model(const ModelOptions& options,
          const Material* mat) const
{
  return NULL;
}


PhysicalModel*
DSSCGeneration::create_boundary_model(const ModelOptions& options,
    const Material* material_A, const Material* material_B) const
{
  return NULL;
}


void
DSSCGeneration::do_setup_solution_variables(void)
{
  declare_solution(Generation, REAL, NODES, "1/(cm^3*s)");
  declare_solution(Distance, REAL, NODES, "cm");
}



void
DSSCGeneration::get_solution_secure(const Elem* elem,
    map<ID, vector<double> >& values, const vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  const NumericVector<Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int d_var = system.variable_number("d");
  //const unsigned int g_var = system.variable_number("G");

  FEType fe_type = system.variable_type(d_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_d;
  //vector<unsigned int> dof_indices_g;

  const vector<vector<Real> >& phi = fe->get_phi();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dof_indices_d, d_var);
  //dof_map.dof_indices(elem, dof_indices_g, g_var);

  unsigned int n_dofs = dof_indices_d.size();

  for (unsigned int n = 0; n < np; n++)
  {
    double d = 0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      d += phi[i][n] * solution(dof_indices_d[i]);
    }

    if (values.count(Distance))
      values[Distance][n] = d;
  }
}



void
DSSCGeneration::do_solve(void)
{

  if (!_d_calculated)
  {
    _calculate_distances();
    _d_calculated = true;
  }



}


void
DSSCGeneration::_calculate_distances(void)
{
  SimulationEnvironment& env = get_environment();
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);
  unsigned short var = system.variable_number("d");

  NumericVector<double>& solution = system.get_solution_vector();

  MeshBase::const_node_iterator it(get_mesh().local_nodes_begin());
  const MeshBase::const_node_iterator end(get_mesh().local_nodes_end());

  for ( ; it != end; ++it)
  {
    const Node* node = *it;
    unsigned int dof = node->dof_number(system.number(), var, 0);

    double d = -1.0;

    SimulationEnvironment::BoundarySideIterator bit(env.boundary_sides_begin());
    const SimulationEnvironment::BoundarySideIterator bend(env.boundary_sides_end());
    for ( ; bit != bend; ++bit)
    {
      const Elem* elem = (bit->first).elem();
      unsigned int s = (bit->first).side();
      if (get_mesh().mesh_dimension() == 1)
      {
        Point dp = elem->point(s) - *node;
        double dd = dp * _direction;
        if (dd >= 0) d = dd;
      }
      else if (get_mesh().mesh_dimension() == 2)
      {

      }
      else // dim = 3
      {

      }

    }

    solution.set(dof, d);

  }

}

