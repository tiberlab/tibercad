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
  // a sensible default
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

  get_parameter("light_direction", _direction);

  // assure it is orthogonal to the "missing" dimensions
  switch (dim)
  {
    case 1:
      _direction(1) = 0;

    case 2:
      _direction(2) = 0;

    default:
      break;
  }

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

  const unsigned int dim = get_mesh().mesh_dimension();

  MeshBase::const_node_iterator it(get_mesh().local_nodes_begin());
  const MeshBase::const_node_iterator end(get_mesh().local_nodes_end());

  for ( ; it != end; ++it)
  {
    const Node& node = *(*it);
    unsigned int dof = node.dof_number(system.number(), var, 0);

    double d = -1.0;

    SimulationEnvironment::BoundarySideIterator bit(env.boundary_sides_begin());
    const SimulationEnvironment::BoundarySideIterator bend(env.boundary_sides_end());
    for ( ; bit != bend; ++bit)
    {
      const Elem* sideelem = (bit->first).elem();
      unsigned int s = (bit->first).side();

      if (dim == 1)
      {
        Point dp = node - sideelem->point(s);
        double dd = dp * _direction;
        if (dd >= 0) d = dd;
      }
      else if (dim == 2)
      {
        AutoPtr<Elem> side_el(sideelem->build_side(s));

        /*
         * have to solve: p0 = node; p1, p2 nodes of side elem
         *
         * -vx * t1 + x0 = x1 + (x2 - x1) * t2
         * -vy * t1 + y0 = y1 + (y2 - y1) * t2
         *
         * t1 >= 0, 0 <= t2 <= 1
         */

        double vx = _direction(0);
        double vy = _direction(1);

        double xa = side_el->point(1)(0) - side_el->point(0)(0);
        double ya = side_el->point(1)(1) - side_el->point(0)(1);

        double xb = node(0) - side_el->point(0)(0);
        double yb = node(1) - side_el->point(0)(1);

        double det = (vy*xa - vx*ya);
        if (det != 0.0)
        {
          double t1 = (xa*yb - xb*ya)/(vy*xa - vx*ya);
          double t2 = (vy*xb - vx*yb)/(vy*xa - vx*ya);

          if ((t1 >= 0) && (t2 >= 0) && (t2 <= 1))
            d = t1;
        }
      }
      else // dim = 3
      {
        AutoPtr<Elem> side_el(sideelem->build_side(s));

        /*
         * have to solve: p0 = node; p1, p2, p3 nodes of side elem (TRI3)
         *
         * -vx * t1 + x0 = x1 + (x2 - x1) * t2 + (x3 - x1) * t3
         * -vy * t1 + y0 = y1 + (y2 - y1) * t2 + (y3 - y1) * t3
         * -vz * t1 + z0 = y1 + (z2 - z1) * t2 + (z3 - z1) * t3
         *
         * t1 > 0, and (p0 - v * t1) inside the side element
         */

        double vx = _direction(0);
        double vy = _direction(1);
        double vz = _direction(2);

        double xa = side_el->point(1)(0) - side_el->point(0)(0);
        double ya = side_el->point(1)(1) - side_el->point(0)(1);
        double za = side_el->point(1)(2) - side_el->point(0)(2);

        double xb = side_el->point(2)(0) - side_el->point(0)(0);
        double yb = side_el->point(2)(1) - side_el->point(0)(1);
        double zb = side_el->point(2)(2) - side_el->point(0)(2);

        double xc = node(0) - side_el->point(0)(0);
        double yc = node(1) - side_el->point(0)(1);
        double zc = node(2) - side_el->point(0)(2);

        double det = (vz*xa*yb - vz*xb*ya - vy*xa*zb + vy*xb*za + vx*ya*zb - vx*yb*za);
        if (det != 0.0)
        {
          double t1 = (xa*yb*zc - xa*yc*zb - xb*ya*zc + xb*yc*za + xc*ya*zb - xc*yb*za) / det;
          double t2 = -(vz*xb*yc - vz*xc*yb - vy*xb*zc + vy*xc*zb + vx*yb*zc - vx*yc*zb) / det;
          double t3 = (vz*xa*yc - vz*xc*ya - vy*xa*zc + vy*xc*za + vx*ya*zc - vx*yc*za) / det;

          if ((t1 >= 0) && (t2 >= 0) && (t3 >= 0))
          {
            Point p(node - t1 * _direction);
            if (side_el->contains_point(p))
              d = t1;
          }

        }


      }

    }

    solution.set(dof, d);

  }

}

