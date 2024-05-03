// $Id$

#include "Poisson.h"
#include "PoissonModel.h"
#include "PoissonBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/plane.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"

// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;
using namespace libMesh;

namespace
{
  void whitney_interpolation(const Elem& elem,
                             const Point& p,
                             vector<double>& w0,
                             vector<Point>& w1)
  {
    unsigned int nn = elem.n_nodes();

    w0.resize(nn);
    w1.resize(elem.n_edges());

    unsigned int dim = 2;

    unique_ptr<FEBase> fe = libMesh::FEBase::build(2, FEType(1, LAGRANGE));
    const vector<vector<Real> >& phi = fe->get_phi();
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();

    vector<Point> pts(1, p);
    fe->reinit(&elem, &pts);

    for (unsigned int i = 0; i < nn; ++i)
      w0[i] = phi[i][0];


    for (unsigned int i = 0; i < nn; ++i)
    {
      unsigned int j = (i + 1) % nn;
      w1[i] = w0[i] * dphi[j][0] - w0[j] * dphi[i][0];
    }
  }
}



Poisson::Poisson(const ModelOptions& options) :
  SimulationInterface(options),
  _my_assembly(this)
{
  // there's nothing to be done
}


Poisson::~Poisson(void)
{
  // there's nothing to be done
}


Poisson*
Poisson::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new Poisson(options);
}



void
Poisson::do_init(void)
{
  parse_options();

  // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  // here we use monomial of constant order, so as to get a
  // finite volume scheme of lowest order
  system.add_variable("u", FIRST);
  system.attach_assemble_object(_my_assembly);
  system.init();
}


void
Poisson::parse_options(void)
{
 myopts.default_boundary_conditions =
   get_option("default_boundary_condition", "zero_field");
 
 string qrule = "trapez";
 qrule = get_option("quadrature_rule", qrule);
 if (qrule == "gauss")
   myopts.quadrature_type = QGAUSS;
 else if (qrule == "trapez")
   myopts.quadrature_type = QTRAP;
 else
   throw InitFailedException("Unknown quadrature rule");

 int intorder = get_option("integration_order", 2);

 if (intorder>7) 
    throw InitFailedException("Invalid integration order");

 myopts.integration_order = static_cast<libMeshEnums::Order>(intorder);
 
 string dualconstr = "barycentric";
 dualconstr = get_option("dual_construction", dualconstr);
 if (dualconstr == "barycentric")
   myopts.dual_constr = BARYCENTER;
 else if (dualconstr == "voronoi")
   myopts.dual_constr = VORONOI;
else
  throw InitFailedException("Unknwon dual construction: " + dualconstr);

string hodgetype = "geometric";
hodgetype = get_option("hodge_type", hodgetype);
if (hodgetype == "geometric")
  myopts.hodge_type = GEOMETRIC;
else if (hodgetype == "Whitney")
  myopts.hodge_type = WHITNEY;
else
  throw InitFailedException("Unknwon Hodge type: " + hodgetype);
}


void
Poisson::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(Potential, REAL, NODES, "V");
  declare_solution(ElField, VECTOR, CELL, "V/cm");
  declare_solution(Displacement, VECTOR, CELL, "C/cm^2");
  declare_solution(Polarization, VECTOR, CELL, "C/cm^2");
  declare_solution(ChargeDensity, REAL, NODES, "C/cm^3");

  // we can define aliases (but the association name -> id
  // has to be surjective)
  add_alias("ElectricField", ElField);
  add_alias("ElPotential", Potential);
}


void
Poisson::do_solve(void)
{

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();


  system.set_options(get_solver_options());
  system.solve();
}


void
Poisson::do_print_info(void)
{
  Messages::info("Module implementation example: Poisson with finite volumes");
}


PhysicalModel*
Poisson::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  return PoissonModel::create(mat, options);
}



PhysicalModel*
Poisson::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return PoissonBoundaryModel::create(boundary, options);
}



void
Poisson::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("u");

  FEType fe_type = system.variable_type(u_var);
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dof_indices, u_var);
  const unsigned int n_dofs = dof_indices.size();

  // cell data variable
  RealGradient field(0);
  RealGradient displacement(0);
  RealGradient pol(0);
  PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

  bool do_displacement = values.count(Displacement);
  bool calculate = false;
  if (do_displacement || values.count(ChargeDensity) || values.count(Polarization))
    calculate = true;

  if (calculate)
    mod.set_element(elem);

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    RealGradient grad(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices[i]);
      grad -= dphi[i][n] * solution(dof_indices[i]);
    }

    field += grad;

    if (values.count(Potential))
      values[Potential][n] = u;


    if (calculate)
    {
      mod.set_point(real_pts[n]);

      mod.calculate();

      const RealTensor& eps = mod.get_permittivity();

      double rho = mod.get_charge_density();

      if (do_displacement)
        displacement -= eps * grad;

      if (values.count(ChargeDensity))
        values[ChargeDensity][n] = rho;
    }
  }

  // this is a very primitive way of estimating the mean value. Better
  // would be to really integrate
  if (do_displacement)
  {
    values[Displacement][0] = displacement(0) / np;
    values[Displacement][1] = displacement(1) / np;
    values[Displacement][2] = displacement(2) / np;
  }

  if ( values.count(Polarization))
  {
    const RealVectorValue& pol = mod.get_polarization();
    values[Polarization][0] = pol(0) / np;
    values[Polarization][1] = pol(1) / np;
    values[Polarization][2] = pol(2) / np;
  }

  if (values.count(ElField))
  {
    values[ElField][0] = field(0) / np;
    values[ElField][1] = field(1) / np;
    values[ElField][2] = field(2) / np;
  }
}



libMesh::Point
Poisson::circumcenter(const libMesh::Elem* elem, int s)
{
  Point x_i(0.0);

  unsigned int dim = elem->dim();

  // (It seems the centroid works better for quadrangles)
  //if (dim == 2)
  if ((dim == 2) && (elem->n_nodes() == 3))
  {
    Point a, b, c;
    
    /*
    if ((s >= 0) && (elem->n_nodes() > 3))
    {
      auto side = elem->side_ptr(s);
      a = side->point(0);
      b = side->point(1);
      Point v1(b - a);

      // look for adjacent side that makes smallest angle
      unsigned int ns = elem->n_sides();

      auto s2 = elem->side_ptr((s+1)%ns);
      Point v2(s2->point(1) - s2->point(0));

      auto s3 = elem->side_ptr((s-1)%ns);
      Point v3(s3->point(0) - s3->point(1));

      double cosa = -(v1 * v2) / (v1.norm() * v2.norm());
      double cosb =  (v1 * v3) / (v1.norm() * v3.norm());

      if (cosa > cosb)
        c = b + v2;
      else
        c = a + v3;
      
    }*/
    //else
    {
      a = elem->point(0);
      b = elem->point(1);
      c = elem->point(2);
    }

    double d = 2 * (a(0) * (b(1) - c(1)) +
                    b(0) * (c(1) - a(1)) + c(0) * (a(1) - b(1)));

    x_i(0) = a.norm_sq() * (b(1) - c(1)) + b.norm_sq() * (c(1) - a(1)) +
             c.norm_sq() * (a(1) - b(1));
    x_i(1) = a.norm_sq() * (b(0) - c(0)) + b.norm_sq() * (c(0) - a(0)) +
             c.norm_sq() * (a(0) - b(0));
    x_i(1) *= -1;
    x_i /= d;

    if (!elem->contains_point(x_i))
    {
      cerr << "outside!\n"; 
      //x_i = elem->centroid();
    }
    else
    {
      // TODO we might check also distance from sides,
      // and change point is to near to an edge
    }
  }
  else if ((dim == 3) && (elem->n_nodes() == 4))
  {
    // tetrahedron
    // circumcenter
    Point u1(elem->point(1) - elem->point(0));
    Point u2(elem->point(2) - elem->point(0));
    Point u3(elem->point(3) - elem->point(0));

    double l1 = u1.norm_sq();
    double l2 = u2.norm_sq();
    double l3 = u3.norm_sq();

    x_i = u2.cross(u3);
    double den = 2 * u1 * x_i;

    x_i *= l1;
    x_i += l2 * u3.cross(u1) + l3 * u1.cross(u2);

    x_i /= den;

    x_i += elem->point(0);

    if (!elem->contains_point(x_i))
    {
      x_i = elem->centroid();
    }
  }
  else
    x_i = elem->centroid();

  return(x_i);
}


void
Poisson::hodge_pd(const libMesh::Elem* elem,
                  libMesh::Point& center,
                  libMesh::DenseMatrix<double>& hodge,
                  vector<pair<unsigned int, unsigned int>>& incidence,
                  vector<double>& vol)
{
  unsigned int dim = elem->dim();
  unsigned int ne = elem->n_edges();
  unsigned int nn = elem->n_nodes();

  hodge.resize(ne, ne);
  incidence.resize(0);
  incidence.reserve(ne);

  vol = vector<double>(nn, 0.0);

  // circumcenter and thus Voronoi-construction works only
  // for triangles, otherwise we fall back to barycentric Hodge
  if ((nn != 3) || (myopts.dual_constr == BARYCENTER))
    center = elem->centroid();
  else
    center = circumcenter(elem);

  vector<Point> prim(ne);

  // get primal 1-chains
  // Works ONLY for 1D and 2D in this form
  for (unsigned int i = 0; i < nn; ++i)
  {
    unsigned int j = (i + 1) % nn;
    prim[i] = (elem->point(j) - elem->point(i));

    incidence.push_back(make_pair(i, j));

    vol[0] = vol[1] = 0.5 * elem->volume();
  }

  if (dim == 1)
  {
    // in this case the Hodge is automatically diagonal
    // dual form has volume 1
    hodge(0, 0) = 1.0 / prim[0].norm();

  }
  else if (dim == 2)
  {

    vector<Point> dual(ne);
    vector<Point> x_m(ne);
    for (unsigned int i = 0; i < nn; ++i)
    {
      unsigned int j = (i + 1) % nn;
      x_m[i] = 0.5 * (elem->point(j) + elem->point(i));

      dual[i] = center - x_m[i];
    }

    auto cross2d = [] (Point& a, Point& b) { return (a(0)*b(1) - a(1)*b(0)); };

    if (myopts.hodge_type == WHITNEY)
    {
      vector<double> w0;
      vector<Point> w1;

      Point p = libMesh::FEInterface::inverse_map(2, libMesh::FEType(), elem, center);
      whitney_interpolation(*elem, p, w0, w1);

      for (unsigned int i = 0; i < nn; ++i)
      {
        hodge(i, 0) += 0.5 * cross2d(w1[0], dual[i]);
        hodge(i, 1) += 0.5 * cross2d(w1[1], dual[i]);
        hodge(i, 2) += 0.5 * cross2d(w1[2], dual[i]);

        double e_x_ed = cross2d(prim[i], dual[i]);

        // calculate the volume contribution associated to the node
        unsigned int r = (i + 1) % nn;
        vol[i] += 0.25 * e_x_ed;
        vol[r] += 0.25 * e_x_ed;
      }

      for (unsigned int i = 0; i < nn; ++i)
      {
        p = libMesh::FEInterface::inverse_map(2, libMesh::FEType(), elem, x_m[i]);
        whitney_interpolation(*elem, p, w0, w1);
        
        hodge(i, 0) += 0.5 * cross2d(w1[0], dual[i]);
        hodge(i, 1) += 0.5 * cross2d(w1[1], dual[i]);
        hodge(i, 2) += 0.5 * cross2d(w1[2], dual[i]);

      }
    }
    else
    {

      for (unsigned int i = 0; i < ne; ++i)
      {
        double e_norm = prim[i].norm_sq();
        double e_x_ed = cross2d(prim[i], dual[i]);
        double e_dot_ed = prim[i] * dual[i];

        // calculate the volume contribution associated to the node
        unsigned int r = (i + 1) % nn;
        vol[i] += 0.25 * e_x_ed;
        vol[r] += 0.25 * e_x_ed;

        hodge(i, i) = e_x_ed / e_norm;

        // the following takes for side i the next two
        // sides for expanding skew-transformed 1-chain:
        // -J e_i = a_i^(i+1)*e_(i+1) + a_i^(i+2)*e_(i+2)
        // It should work for all n-angles, unless two
        // subsequent sides are collinear
        for (unsigned int j = 1; j < 3; ++j)
        {
          unsigned int k = (i + j) % ne;
          unsigned int l = (i + 3 - j) % ne;

          double ei_dot_el = prim[i] * prim[l];
          double ek_x_el = cross2d(prim[k], prim[l]);

          hodge(i, k) = (e_dot_ed / e_norm) * (ei_dot_el / ek_x_el);
        }
      }
    }
  }
}



void
Poisson::assemble(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  
  get_scaling().set_length_scaling(get_mesh_units());
  double x0 = get_mesh_units(); 
  double Lambda = Constants::e * 1e6 * x0*x0 / Constants::e0;
  
  DofMap& dof_map =  system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  
  // The Hodge dual
  DenseMatrix<Number> H;
  // The discrete exterior derivative D0
  DenseMatrix<Number> D0;


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    // number of neighbors (= number of sides, usually)
    unsigned int n_neigh = elem->n_neighbors();

    // all DoF indices in the current block
    // in this case DoFs on all nodes
    vector<unsigned int> dof_indices;

    dof_map.dof_indices(elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    D0.resize(n_dofs, n_dofs);

    Point center;
    vector<double> vol;
    
    // the incidence contains the pairs of indices of nodes
    // in the same order as the entries in H
    vector<pair<unsigned int, unsigned int>> inc;
    hodge_pd(elem, center, H, inc, vol);


    //cerr << "H = \n";
    //cerr << H << endl;
    //for (unsigned int i = 0; i < vol.size(); ++i)
    //  cerr << vol[i] / elem->volume() << " ";
    //cerr << endl << endl;


    PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

    mod.set_element(elem);
    mod.set_point(center);
    mod.calculate();

    // permittivity
    const RealTensor &eps_i = mod.get_permittivity();
    // for now assume eps_k = e_i * I
    double e_i = eps_i.tr() / 3.0;

    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
    {
      unsigned int j = (i + 1) % elem->n_nodes();
      D0(i, i) = -1;
      D0(i, j) =  1;

      mod.set_point(elem->point(i));
      mod.calculate();

      // polarization
      const RealVectorValue &pol_i = mod.get_polarization();

      double rho = mod.get_charge_density() * Lambda;
      Fe(i) += rho * vol[i];

    }

    H.right_multiply(D0);
    //cerr << "H*D = \n";
    //cerr << H << "\n";
    H *= -e_i;

    // loop over the neighbors
    // j is used as column index in the Ke matrix
    for (unsigned int k = 0, j = 0; k < n_neigh; ++k)
    {

      // interface model if present
      PoissonBoundaryModel *mod_int =
          get_interface_model<PoissonBoundaryModel>(elem, k);

      if (mod_int != nullptr)
      {
        auto side_elem = elem->side_ptr(k);
        double side_vol = side_elem->volume();
        Point side_ctr = side_elem->centroid();

        mod_int->calculate(elem, k, side_ctr);

        double a, b, c;
        mod_int->get_coefficients(a, b, c);

        // we use a penalty approach here for its simplicity
        if ((b < 1e-10) && (b >= 0))
          b = 1e-10;
        else if ((b > -1e-10) && (b <= 0))
          b = -1e-10;

        a /= b;
        c /= b;

        for (unsigned int i = 0; i < elem->n_nodes(); ++i)
        {
          if (elem->is_node_on_side(i, k))
          {
            Ke(i, i) += a * side_vol;
            Fe(i) += c * side_vol;
          }
        }
      }
    }

    for (int i = 0; i < H.n(); ++i)
    {
      int ii = (i - 1) % (int) H.n();
      if (ii < 0) ii += H.n();

      for (unsigned int j = 0; j < H.n(); ++j)
      {
        Ke(i, j) += H( i, j);
        Ke(i, j) -= H(ii, j);
      }
    }

    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);


  }
  system.matrix->close();
  //system.matrix->print_matlab("K.m");
  system.rhs->close();
  //system.rhs->print_matlab("F.m");

}
