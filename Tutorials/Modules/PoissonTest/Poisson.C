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

// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;
using namespace libMesh;




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
  system.add_variable("u", CONSTANT, MONOMIAL);
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
  UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));

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



void
Poisson::assemble(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  
  get_scaling().set_length_scaling(get_mesh_units());
  double x0 = get_mesh_units(); 
  double Lambda = Constants::e * 1e6 * x0*x0;
  
  DofMap& dof_map =  system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    // number of neighbors (= number of sides, usually)
    unsigned int n_neigh = elem->n_neighbors();

    // all DoF indices in the current block
    vector<unsigned int> dof_indices;

    // only DoFs of the current element
    vector<unsigned int> rows;

    dof_map.dof_indices(elem, rows);
    // NOTE: as we hardcoded 0th order monomials rows will
    // always be of size 1
    dof_indices.push_back(rows[0]);

    for (unsigned int i = 0; i < n_neigh; ++i)
    {
      const Elem* neigh = elem->neighbor_ptr(i);

      vector<unsigned int> tmpdofs;

      if (neigh != nullptr)
      {
        dof_map.dof_indices(neigh, tmpdofs);
        dof_indices.push_back(tmpdofs[0]);
      }
    }

    const unsigned int n_dofs = dof_indices.size();

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(1, n_dofs);
    Fe.resize(1);

    PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

    mod.set_element(elem);

    // calculate parameters we need from this element
    Point x_i(elem->centroid());

    double vol_i = elem->volume();

    mod.set_point(x_i);
    mod.calculate();

    // permittivity
    const RealTensor& eps_i = mod.get_permittivity()*Constants::e0;
    // for now assume eps_k = e_i * I
    double e_i = eps_i.tr() / 3.0;

    // polarization
    const RealVectorValue& pol_i = mod.get_polarization();

    double rho =  mod.get_charge_density() * Lambda;
    Fe(0) += rho * vol_i;

    // loop over the neighbors
    // j is used as column index in the Ke matrix
    for (unsigned int k = 0, j = 0; k < n_neigh; ++k)
    {
      const Elem* neigh = elem->neighbor_ptr(k);

      // interface model if present
      PoissonBoundaryModel* mod_int =
        get_interface_model<PoissonBoundaryModel>(elem, k);

      // First set up some geometric quantities.
      //
      // We need the distance from the element interface
      // (from both sides), and the area of the side,
      // and the projection on the side.
      // In the following m indicates the middle = interface
      Point x_m;
      Point x_k;
      double h_im, h_mk = 0.0;
      double A_ik = 1;

      // the normal, pointing out of the current element (i)
      Point n_ik;


      if (neigh != nullptr)
      {
         // advance index counter
         ++j;

         x_k = neigh->centroid();

      }
      else if (mod_int == nullptr)
      {
        // in this case there is nothing to be done
        // (natural boundary conditions)
        continue;
      }

      if (dim == 1)
      {
        x_m = elem->point(k);
        Point d(x_m);
        d -= x_i;
        h_im = d.norm();

        n_ik = d / h_im;

        if (neigh != nullptr)
        {
          Point d(x_k);
          d -= x_m;
          h_mk = d.norm();
        }
      }

      if (dim == 2)
      {
        // get the projection onto a side
        // we need the side
        // NOTE: should be adapted for any orientation in 3D
        auto side = elem->side_ptr(k);
        Point x1(side->point(0));
        Point x2(side->point(1));

        Point d(x2 - x1);
        Point y = x1 - x_i;

        double det = -(d(0)*d(0) + d(1)*d(1));
        double t = (d(1)*y(0) - d(0)*y(1)) / det;

        n_ik(0) =  d(1);
        n_ik(1) = -d(0);
        n_ik *= t;
        h_im = n_ik.norm();
        n_ik = n_ik / h_im;

        A_ik = side->volume();

        if (neigh != nullptr)
        {
           y = x1 - x_k;

           double det = -(d(0) * d(0) + d(1) * d(1));
           double t = (d(1) * y(0) - d(0) * y(1)) / det;

           Point n_ki;
           n_ki(0) = d(1);
           n_ki(1) = -d(0);
           n_ki *= t;
           h_mk = n_ki.norm();
        }
      }

      if (dim == 3)
      {
        // get the projection onto a side
        // we need the side
        auto side = elem->side_ptr(k);
        Point x1(side->point(0));
        Point x2(side->point(1));
        Point x3(side->point(2));

        libMesh::Plane p(x1, x2, x3);
        Point pp = p.closest_point(x_i);

        n_ik = pp - x_i;
        
        h_im = n_ik.norm();
        n_ik = n_ik / h_im;

        A_ik = side->volume();

        if (neigh != nullptr)
        {
          pp = p.closest_point(x_k);

          Point n_ki = pp - x_k;
        
          h_mk = n_ki.norm();
        }
      }


      double a = 0.0, b = 0.0, c = 0.0;
      if (mod_int != nullptr)
      {
        mod_int->calculate(elem, k, x_m);
        mod_int->get_coefficients(a, b, c);
      }

      // phi_m is linear function in phi_i and phi_k:
      // phi_m = c_i * phi_i + c_k * phi_k + f
      double c_i = 0;
      double c_k = 0;
      double f = 0;


      // interface charge, if present
      double sigma_int = 0.0;


      if ((b == 0.0) && (a != 0.0))
      {
        // this is a Dirichlet BC
        f = c / a;
      }
      else if (neigh != nullptr)
      {

        PoissonModel &mod = *get_bulk_model<PoissonModel>(neigh);

        mod.set_element(neigh);

        mod.set_point(x_k);
        mod.calculate();

        const RealTensor &eps_k = mod.get_permittivity() * Constants::e0;
        // for now assume eps_k = e_k * I
        double e_k = eps_k.tr() / 3.0;
        const RealVectorValue &pol_k = mod.get_polarization();

        // polarization induced interface charge density
        double sigma_pol = (pol_k - pol_i) * n_ik;

        // interface potential discontinuity (due to e.g. dipole)
        // double delta_phi = 0.0;

        double denom = h_mk * e_i + h_im * e_k;
        c_i = h_mk * e_i / denom;
        c_k = h_im * e_k / denom;
        f = /* -h_im*e_k / denom * delta_phi */
            -h_im * h_mk / denom * (sigma_pol + sigma_int);

      }


      // now we use
      //   D_i = -e_i * (phi_m - phi_i) / h_im + P_i*n_ik - 0.5*sigma_int;
      Ke(0, 0) += e_i / h_im * (1 - c_i) * A_ik;

      if (neigh != nullptr)
        Ke(0, j) -= e_i / h_im * c_k * A_ik;

      Fe(0) -= (pol_i * n_ik - 0.5 * sigma_int - e_i / h_im * f) * A_ik;

    }

    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, rows, dof_indices);
    system.rhs->add_vector(Fe, rows);

  }
  //system.matrix->close();
  //system.matrix->print_matlab("K.m");
  //system.rhs->close();
  //system.rhs->print_matlab("F.m");

}
