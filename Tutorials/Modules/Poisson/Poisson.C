// $Id$

#include "Poisson.h"
#include "PoissonModel.h"
#include "PoissonBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"

#include "equation_systems.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"


// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;


Poisson*
Poisson::_this = NULL;


Poisson::Poisson(const ModelOptions& options) :
  SimulationInterface(options)
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
  // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  system.add_variable("u", FIRST);
  system.attach_assemble_function(assemble);
  system.init();
}


void
Poisson::parse_options(void)
{
 myopts.default_boundary_conditions =
   get_options().get_option("default_boundary_condition", "zero_field");
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
}


void
Poisson::do_solve(void)
{
  // this is dirty, but at the moment we have to provide a static method
  // for assembly. Thsi should change in future
  _this = this;

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();


  system.set_options(get_solver_options());
  system.solve();
}


void
Poisson::do_print_info(void)
{
  Messages::info("Module implementation example: simple Poisson");
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

  const NumericVector<Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("u");

  FEType fe_type = system.variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

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
      mod.set_point(p[n]);

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
Poisson::do_assemble(EquationSystems& es, const std::string& system_name)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  const unsigned int uvar = system.variable_number("u");

  FEType fe_type = dof_map.variable_type(uvar);

  // the volume finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, SECOND);
  fe->attach_quadrature_rule(&qrule);

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // the surface finite element
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qface(dim - 1, THIRD);
  fe_face->attach_quadrature_rule(&qface);

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();

  vector<unsigned int> dof_indices;

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

    mod.set_element(elem);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      mod.set_point(q_point[qp]);

      mod.calculate();

      const RealTensor& eps = mod.get_permittivity();
      const RealVectorValue& pol = mod.get_polarization();
      double rho = Constants::e * mod.get_charge_density();



      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
          Ke(i, j) += JxW[qp] * Constants::e0 * (dphi[i][qp] * (eps * dphi[j][qp]));

	Fe(i) += JxW[qp] * (rho * phi[i][qp] + pol * dphi[i][qp]);
      }

    }



    // the sides
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      PoissonBoundaryModel* mod_int =
        get_interface_model<PoissonBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        fe_face->reinit(elem, s);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod_int->calculate(elem, s, qface_point[qp]);

          double a, b, c;
          mod_int->get_coefficients(a, b, c);

          // we use a penalty approach here for its simplicity
          if ((b < 1e-10) && (b >= 0)) b = 1e-10;
          else if ((b > -1e-10) && (b<= 0)) b = -1e-10;

          a /= b;
          c /= b;

          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
              Ke(i, j) += a * JxW_face[qp] * (phi_face[i][qp] * phi_face[j][qp]);

            Fe(i) += c * JxW_face[qp] * phi_face[i][qp];
          }
        }
      }
    }

    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);

  }
  //system.matrix->close();
  //system.matrix->print_matlab("K.m");
  //system.rhs->close();
  //system.rhs->print_matlab("F.m");

}
