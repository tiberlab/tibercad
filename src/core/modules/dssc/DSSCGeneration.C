// $Id$

#include "DSSCGeneration.h"
#include "TiberLinearSystem.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Database.h"

#include "TiberModule.h"

#include "mesh_base.h"
#include "dof_map.h"

#include <limits>


using namespace std;


DSSCGeneration::DSSCGeneration(const ModelOptions& options) :
  SimulationInterface(options),
  _d_calculated(false),
  _direction(0),
  _intensity(0),
  _alpha(-1)
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

  // light intensity is given in x sun
  get_parameter("light_intensity", _intensity);

  if (get_options().find_option("dye"))
    _read_spectrum();
  else
  {
    // alpha is given in um^-1
    _alpha = 0.15;
    get_parameter("alpha", _alpha);
  }


  // length units in cm
  double mesh_units = 100 * get_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  //system.attach_assembly_routine(assemble_system);


  system.add_variable("d", libMeshEnums::FIRST);

  //system.add_vector("G");

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

    if (d < 0) { d = -d; }

    if (values.count(Distance))
      values[Distance][n] = d;

    if (values.count(Generation))
    {
      double generation = 0.0;
      if (_alpha >= 0)
       generation = 1e4 * 1.5e17 * _intensity * _alpha * exp(-d * 1e4 *_alpha);
//       generation = 1e4 * 1.0e17 * _intensity * _alpha * exp(-d * 1e4 *_alpha);
//       generation = 1e4 * 1.5e17 * _intensity * _alpha * exp(-d * _alpha);
      else
      {
        generation = _intensity * _integrate(d);
      }

      values[Generation][n] = generation;
    }
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
  unsigned short sysnr = system.number();
  unsigned short var = system.variable_number("d");

  NumericVector<double>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const double x0 = get_scaling().get_calc_mesh_units();

  MeshBase::const_node_iterator it(get_mesh().local_nodes_begin());
  const MeshBase::const_node_iterator end(get_mesh().local_nodes_end());

  for ( ; it != end; ++it)
  {
    const Node& node = *(*it);
    unsigned short ncomp = node.n_comp(sysnr, var);
    if (ncomp == 0) continue;

    unsigned int dof = node.dof_number(sysnr, var, 0);

    double d = numeric_limits<double>::max();

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
        d = dd < d ? dd : d;
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
            d = t1 < d ? t1 : d;
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

          //if ((t1 >= 0) && (t2 >= 0) && (t3 >= 0))
          if (t1 >= 0)
          {
            Point p(node - t1 * _direction);
            if (side_el->contains_point(p))
              d = t1 < d ? t1 : d;
          }

        }


      }

    }

    solution.set(dof, x0 * d);

  }

}


void
DSSCGeneration::_read_spectrum(void)
{
  // we read the solar spectrum and the absorption spectrum, and
  // then we interpolate the latter at the wavelength given in
  // the former.


  // a line buffer
  const size_t buf_len = 256;
  char buf[buf_len];

  Database db(get_option("dye", ""));
  ifstream is(db.get_data_file().c_str());
  if (is.fail() || !is.good())
    throw InitFailedException("Cannot read absorption "
        "spectrum from file " + db.get_data_file());

  double gap = 0;

  size_t i = 0;
  vector<double> l_tmp;
  vector<double> a_tmp;
  while (is.good())
  {
    if (i == l_tmp.size())
    {
      size_t n_new = l_tmp.size() + 100;
      l_tmp.reserve(n_new);
      a_tmp.reserve(n_new);
    }

    is.getline(buf, buf_len);
    if (buf[0] != '#')
    {
      istringstream in(buf);

      double l, s;
      if (in >> l)
      {
        if (in >> s)
        {
          l_tmp.push_back(l);
          // alpha given in um^-1
          a_tmp.push_back(1e4 * s);
          i++;
        }
        else
        {
          gap = Constants::e * l;
        }
      }
    }
  }
  is.close();

  l_tmp.resize(l_tmp.size());
  a_tmp.resize(a_tmp.size());

  if (gap == 0)
    throw InitFailedException("No energy gap specified "
        "in file " + db.get_data_file());


  // now read the solar spectrum

  db.set_material("Sun1p5am", get_option("illumination_spectrum", ""));
  is.open(db.get_data_file().c_str());

  if (is.fail() || !is.good())
    throw InitFailedException("Cannot read spectrum "
        "from file " + db.get_data_file());

  i = 0;
  while (is.good())
  {
    if (i == _lambda.size())
    {
      size_t n_new = _lambda.size() + 100;
      _lambda.reserve(n_new);
      _spectrum.reserve(n_new);
    }

    is.getline(buf, buf_len);
    if (buf[0] != '#')
    {
      istringstream in(buf);

      double l, s;
      if (in >> l >> s)
      {
        _lambda.push_back(l);
        _spectrum.push_back(1e-7 * s / gap);
        i++;
      }
    }
  }
  is.close();

  _lambda.resize(_lambda.size());
  _spectrum.resize(_spectrum.size());



  _absorption.resize(_lambda.size());
  size_t j = 0;
  size_t jmin = 0;
  size_t k = 0;
  while (j < _lambda.size())
  {
    while ((k < l_tmp.size()) && (l_tmp[k] < _lambda[j]))
      k++;

    if (k == l_tmp.size()) break;

    if (k == 0)
    {
      jmin = j;
      _absorption[j] = 0;
    }
    else
    {
      double dl = l_tmp[k] - l_tmp[k - 1];
      double x  = (_lambda[j] - l_tmp[k - 1]) / dl;

      _absorption[j] = x * a_tmp[k] + (1 - x) * a_tmp[k - 1];
    }
    j++;
  }

  // limit to the part that is non zero
  l_tmp.resize(j - jmin);
  a_tmp.resize(j - jmin);
  vector<double> s_tmp(j - jmin);
  for (size_t i = 0; i < l_tmp.size(); i++)
  {
    size_t idx = jmin + i;
    l_tmp[i] = _lambda[idx];
    a_tmp[i] = _absorption[idx];
    s_tmp[i] = _spectrum[idx];
  }

  _lambda = l_tmp;
  _absorption = a_tmp;
  _spectrum = s_tmp;

  //for (size_t i = 0; i < l_tmp.size(); i++)
  //  cerr << _lambda[i] << " " << _spectrum[i] << "  " << _absorption[i] << endl;
}



double
DSSCGeneration::_integrate(double d)
{
  double gen = 0.0;

  size_t N = _lambda.size();

  for (size_t i = 0; i < N - 1; i++)
  {
    double dl = _lambda[i + 1] - _lambda[i];

    double a1 = _absorption[i + 1];
    double f1 = _spectrum[i + 1] * a1 * exp(-a1 * d);

    double a0 = _absorption[i];
    double f0 = _spectrum[i] * a0 * exp(-a0 * d);

    gen += 0.5 * dl * (f1 + f0);
  }

  return gen;
}
