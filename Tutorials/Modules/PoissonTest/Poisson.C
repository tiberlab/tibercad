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
  
  get_scaling().set_length_scaling(1.0);
  double Lambda = Constants::e * 1e6;
  
  DofMap& dof_map =  system.get_dof_map();

  /*
  const unsigned int uvar = system.variable_number("u");

  FEType fe_type = dof_map.variable_type(uvar);

  // the finite element
  UniquePtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  UniquePtr<QBase> qrule(QBase::build(myopts.quadrature_type, dim, myopts.integration_order));
  fe->attach_quadrature_rule(qrule.get());

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  // the surface finite element
  UniquePtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  UniquePtr<QBase> qface(QBase::build(myopts.quadrature_type, dim-1, myopts.integration_order));
  fe_face->attach_quadrature_rule(qface.get());

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();
  */

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    // number of neighbors (= number of sides, usually)
    unsigned int n_neigh = elem->n_neighbors();

    // list of neighbor element pointers
    std::vector<const Elem*> neighbors(n_neigh, nullptr);

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
        neighbors[i] = neigh;
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

    mod.set_point(x_i);
    mod.calculate();

    const RealTensor& eps_i = mod.get_permittivity()*Constants::e0;
    const RealVectorValue& pol_i = mod.get_polarization();

    double rho =  mod.get_charge_density() * Lambda;
    Fe(0) += rho * elem->volume();

    // loop over the neighbors
    for (unsigned int k = 0; k < n_neigh; ++k)
    {
      const Elem* neigh = elem->neighbor_ptr(k);

      if (neigh == nullptr)
        continue;

      PoissonModel &mod = *get_bulk_model<PoissonModel>(neigh);

      mod.set_element(neigh);

      // calculate parameters we need from this element
      Point x_k(neigh->centroid());

      mod.set_point(x_k);
      mod.calculate();

      const RealTensor &eps_k = mod.get_permittivity() * Constants::e0;
      const RealVectorValue &pol_k = mod.get_polarization();

      // now we need the distance from the element interface from both sides,
      // and the area of the side
      // In the following m indicates the middle = interface

      double h_im, h_mk;
      double A_ik = 1;

      if (dim == 1)
      {
        Point d(elem->point(k));
        d -= x_i;
        h_im = d.norm();

        d = x_k;
        d -= elem->point(k);
        h_mk = d.norm();
        
      }

      if (dim == 2)
      {
        // get the projection onto a side
      }

      if (dim == 3)
      {
        // get the projection onto a side
      }

      Ke(0, 0) += 1;

    }

      //PoissonBoundaryModel* mod_int =
       // get_interface_model<PoissonBoundaryModel>(elem, s);

      //if (mod_int != NULL)
 

    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, rows, dof_indices);
    system.rhs->add_vector(Fe, rows);

  }
  system.matrix->close();
  system.matrix->print_matlab("K.m");
  //system.rhs->close();
  //system.rhs->print_matlab("F.m");

}
