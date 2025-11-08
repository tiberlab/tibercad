// $Id$

#include "Poisson.h"
#include "PoissonModel.h"
#include "PoissonBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "DECdual.h"
#include "PolyTope.h"

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

#include <queue>

// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;
using namespace libMesh;

namespace
{
  /*
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
  */
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
  system.add_variable("u", CONSTANT, MONOMIAL);
  system.attach_assemble_object(_my_assembly);
  system.init();
}


void
Poisson::parse_options(void)
{
 myopts.default_boundary_conditions =
   get_option("default_boundary_condition", "zero_field");
 
}


void
Poisson::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(Potential, REAL, CELL, "V");
  declare_solution(ElField, VECTOR, CELL, "V/cm");
  declare_solution(Displacement, VECTOR, CELL, "C/cm^2");
  declare_solution(Polarization, VECTOR, CELL, "C/cm^2");
  declare_solution(ChargeDensity, REAL, CELL, "C/cm^3");

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
  Messages::info("Module implementation example: Poisson with dual DEC");
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
    RealGradient grad(0);

    double u = solution(dof_indices[0]);
    //grad -= dphi[i][n] * solution(dof_indices[i]);

    //field += grad;

    if (values.count(Potential))
      values[Potential][n] = u;

    if (calculate)
    {
      mod.set_point(real_pts[n]);

      mod.calculate();

      const RealTensor &eps = mod.get_permittivity();

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





void
Poisson::hodge_pd(const libMesh::Elem* elem,
                  libMesh::Point& center,
                  libMesh::DenseMatrix<double>& hodge,
                  vector<pair<unsigned int, unsigned int>>& incidence,
                  vector<double>& vol)
{
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

    // number of actual neighbors
    unsigned int n_neigh = 0;
    for (unsigned int i = 0; i < elem->n_sides(); ++i)
    {
      if (elem->neighbor_ptr(i) != nullptr) ++n_neigh;
    }

    // all DoF indices in the current block
    // in this case DoFs on all nodes
    vector<unsigned int> dof_indices;

    dof_map.dof_indices(elem, dof_indices);
    vector<unsigned int> dof_indices_all(dof_indices);
    dof_indices_all.reserve(n_neigh + 1); // prepare capacity

    // this is = 1 in the present case
    const unsigned int n_dofs = dof_indices.size();

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(n_dofs, n_neigh + 1);
    Fe.resize(n_dofs);

    D0.resize(n_neigh, n_neigh + 1);
    H.resize(n_dofs, n_neigh);

    Point center = elem->vertex_average();
    double vol = elem->volume();
    
    PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

    mod.set_element(elem);
    mod.set_point(center);
    mod.calculate();

    // permittivity
    const RealTensor &eps_i = mod.get_permittivity();
    // for now assume eps_k = e_i * I
    double e_i = eps_i.tr() / 3.0;

    // polarization
    const RealVectorValue &pol_i = mod.get_polarization();

    double rho = mod.get_charge_density() * Lambda;
    Fe(0) = rho * vol;

    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
    {
      DECdual dual(elem, i);
      PolyTope pt;
      dual.get_polytope(pt);

    }

    // loop over all sides, but do something only on sides with neighbor
    for (unsigned int i = 0, j = 0; i < elem->n_sides(); ++i)
    {
      const Elem* neigh = elem->neighbor_ptr(i);

      auto side = elem->side_ptr(i);

      // the primal edge
      Point eij = side->point(0) - side->point(1);

      auto cross2d = [] (Point& a, Point& b) { return (a(0)*b(1) - a(1)*b(0)); };
      
      // interface model if present
      PoissonBoundaryModel *mod_int =
          get_interface_model<PoissonBoundaryModel>(elem, i);

      if (mod_int != nullptr)
      {
        Point eij_dual = side->vertex_average();

        mod_int->calculate(elem, i, eij_dual);

        eij_dual -= center;
        double hij2 = eij_dual.norm_sq();

        double value = std::abs(cross2d(eij_dual, eij)) / hij2;


        double a, b, c;
        mod_int->get_coefficients(a, b, c);

        // we use a penalty approach here for its simplicity
        if ((b < 1e-10) && (b >= 0))
          b = 1e-10;
        else if ((b > -1e-10) && (b <= 0))
          b = -1e-10;

        a /= b;
        c /= b;

        Ke(0, 0) += a * value;
        Fe(0) += c * value;
      }

      // if there is no neighbor element, we can skip the rest
      if (neigh == nullptr) continue;

      vector<unsigned int> dof_indices_neigh;

      // discrete derivation
      D0(j, 0) = -1;
      D0(j, j+1) =  1;

      dof_map.dof_indices(neigh, dof_indices_neigh);
      // here we rely on having a single DOF per element
      dof_indices_all.push_back(dof_indices_neigh[0]);

      // now the discrete Hodge

      Point eij_dual = neigh->vertex_average();
      eij_dual -= center;
      double hij2 = eij_dual.norm_sq();

      H(0, j) = std::abs(cross2d(eij_dual, eij)) / hij2;

      ++j;
    }

    H.right_multiply(D0);
    H *= -e_i;
    //cerr << "H*D = \n";
    //cerr << H << "\n";
    Ke += H;


    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices, dof_indices_all);
    system.rhs->add_vector(Fe, dof_indices);


  }
  system.matrix->close();
  system.matrix->print_matlab("K.m");
  system.rhs->close();
  //system.rhs->print_matlab("F.m");

}
