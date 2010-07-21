// $Id$

#include "MyPoisson.h"
#include "MyPoissonModel.h"
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
// The first string is the class name of the object to be created,
// the second one is the name of the module as it should be referred
// in the input file (the Makefile defines MODULE_NAME, which can be used here).
TIBER_MODULE(MyPoisson, MODULE_NAME)


using namespace std;


MyPoisson*
MyPoisson::_this = NULL;


MyPoisson::MyPoisson(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}


MyPoisson::~MyPoisson(void)
{
  // there's nothing to be done
}


MyPoisson*
MyPoisson::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new MyPoisson(options);
}



void
MyPoisson::do_init(void)
{
  TiberLinearSystem* system = TiberLinearSystem::create(get_equation_systems(),
      get_equation_system_name(), get_solver_options());

  system->add_variable("u", FIRST);
  system->attach_assemble_function(assemble);
  system->init();
}


void
MyPoisson::parse_options(void)
{
}


void
MyPoisson::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(Potential, REAL, NODES, "V");
  declare_solution(Field, VECTOR, CELL, "V/cm");
  declare_solution(Displacement, VECTOR, CELL, "C/cm^2");
  declare_solution(ChargeDensity, REAL, NODES, "C/cm^3");

  // we can define aliases (but the association name -> id
  // has to be surjective)
  add_alias("ElectricField", Field);
}


void
MyPoisson::do_solve(void)
{
  _this = this;

  EquationSystems& es = get_equation_systems();

  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>(get_equation_system_name());


  system.set_options(get_solver_options());
  system.solve();
}


void
MyPoisson::do_print_info(void)
{
  Messages::info("Module implementation example");
}


PhysicalModel*
MyPoisson::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  return MyPoissonModel::create(options);
}



PhysicalModel*
MyPoisson::create_boundary_model(const ModelOptions& options,
    const Material* material_A, const Material* material_B) const
{
  return PoissonBoundaryModel::create(options);
}



void
MyPoisson::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem* system;
  system = &get_equation_systems().get_system<TiberLinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("u");

  FEType fe_type = system->variable_type(u_var);
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

  MyPoissonModel& mod = *get_bulk_model<MyPoissonModel>(elem);

  bool do_displacement = values.count(Displacement);
  bool calculate = false;
  if (do_displacement || values.count(ChargeDensity))
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
      const RealVectorValue& pol = mod.get_polarization();
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

  if (values.count(Field))
  {
    values[Field][0] = field(0) / np;
    values[Field][1] = field(1) / np;
    values[Field][2] = field(2) / np;
  }
}



void
MyPoisson::do_assemble(EquationSystems& es, const std::string& system_name)
{
  TiberLinearSystem& system = static_cast<TiberLinearSystem&>(
      get_equation_systems().get_system(get_equation_system_name()));

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  const unsigned int uvar = system.variable_number("u");

  FEType fe_type = dof_map.variable_type(uvar);


  // the volume finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // the surface finite element
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qface(dim - 1, SIXTH);
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

    MyPoissonModel& mod = *get_bulk_model<MyPoissonModel>(elem);

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
      PoissonBoundaryModel* mod =
        get_interface_model<PoissonBoundaryModel>(elem, s);

      if (mod != NULL)
      {
        fe_face->reinit(elem, s);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod->calculate(elem, s, qface_point[qp]);

          double a, b, c;
          mod->get_coefficients(a, b, c);

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
