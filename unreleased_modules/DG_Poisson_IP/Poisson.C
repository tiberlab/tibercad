// Poisson.C Created on: Feb 11, 2016 Author: mpatria

//DGFEM to solve Poisson Equation

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
#include "fe_interface.h"
#include "fe.h"

// This is needed in order to create the shared module library
#include "TiberModule.h"

// C++ includes
#include <fstream>
#include <iostream>


using namespace std;


Poisson*
Poisson::_this = NULL; //not have exact address to be assigned


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
  parse_options();

  ID  dim = get_mesh().mesh_dimension();

  // create a linear equation system
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  system.add_variable("u", FIRST, MONOMIAL);
  system.attach_assemble_function(assemble); //Register a user function to use in assembling the system matrix and RHS.
  system.init();
}


libMesh::NumericVector<double>&
Poisson::do_get_solution_vector(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
  system.get_solution_vector().close();
  return system.get_solution_vector();
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
  declare_solution(Potential, REAL, CELL, "V");
  declare_solution(ElField, VECTOR, CELL, "V/cm");
  declare_solution(Displacement, VECTOR, CELL, "C/cm^2");
  declare_solution(Polarization, VECTOR, CELL, "C/cm^2");
  declare_solution(ChargeDensity, REAL, NODES, "C/cm^3");
}



void
Poisson::do_solve(void)
{
  // this is dirty, but at the moment we have to provide a static method
  // for assembly. Thsi should change in future
  _this = this;

  libMesh::EquationSystems& es = get_equation_systems();

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
Poisson::get_solution_secure(const Elem* elem, std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const libMesh::NumericVector<Number>& solution = system.get_solution_vector();

  const libMesh::DofMap& dof_map = system.get_dof_map();

  ID dim = get_mesh().mesh_dimension();
  double mesh_units = get_mesh_units();

  const unsigned int u_var = system.variable_number("u");

  libMesh::FEType fe_type = system.variable_type(u_var);
  libMesh::AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  libMesh::AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type));
  libMesh::AutoPtr<FEBase> fe_face_neig(build_finite_element(dim, fe_type));

  libMesh::QGauss qface(dim - 1, CONSTANT);
  fe_face->attach_quadrature_rule(&qface);

  vector<unsigned int> dof_indices;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& q_point = fe->get_xyz();

  const vector<vector<Real> >& phi_face = fe_face->get_phi();
  const vector<Point>& normal_face = fe_face->get_normals();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<Real>& JxW_face = fe_face->get_JxW();

  const vector<vector<Real> >& phi_face_neig = fe_face_neig->get_phi();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dof_indices, u_var);
  const unsigned int n_dofs = dof_indices.size();

  unsigned int n_sides = elem->n_sides();

  RealGradient field_cons(0);
  RealGradient displac(0);
  RealGradient lifting_d(0);
  RealGradient lifting_e(0);
  RealGradient pol(0);


  PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);

  bool do_displacement = values.count(Displacement);
  bool calculate = false;
  if (do_displacement || values.count(ChargeDensity) || values.count(Polarization))//|| values.count(Polarization))
    calculate = true;

  if (calculate)
    mod.set_element(elem);

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    RealGradient grad(0);

    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices[i]);
      grad -= dphi[i][n] * solution(dof_indices[i]);
    }

    field_cons += grad;

    if (values.count(Potential))
      values[Potential][n] = u;

    if (calculate)
    {
      mod.set_point(p[n]);

      mod.calculate(elem, p[n]);

      const libMesh::RealTensor& eps_elem = mod.get_permittivity();

      double rho = mod.get_charge_density();

      if (do_displacement)
        displac -= eps_elem * grad;

      if (values.count(ChargeDensity))
        values[ChargeDensity][n] = rho;
    }
  }

  if (do_displacement)
  {
    double my_det_average,
           a11, a12, a13,
           a21, a22, a23,
           a31, a32, a33;

    for (unsigned int s = 0; s < n_sides; s++)
    {
      const Elem* elem_neig = elem->neighbor(s);

      if (elem_neig != NULL)
      {
        fe_face->reinit(elem, s);

        PoissonModel& mod_neig = *get_bulk_model<PoissonModel>(elem_neig);

        std::vector<Point> qface_point_neig;

        AutoPtr<Elem> elem_side(elem->build_side(s));
        const double h_elem = (elem->volume() / elem_side->volume()) * mesh_units;

        unsigned int penalty = 4;

        std::vector<unsigned int> neig_dof_indices;
        dof_map.dof_indices (elem_neig,neig_dof_indices);

        FEInterface::inverse_map(elem->dim(), fe_type, elem_neig, qface_point, qface_point_neig);
        fe_face_neig->reinit(elem_neig, &qface_point_neig);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod.set_point(qp);

          mod.calculate(elem, qface_point[qp]);
          const libMesh::RealTensor& eps_elem = mod.get_permittivity();

          mod_neig.calculate(elem_neig, qface_point_neig[qp]);
          const libMesh::RealTensor& eps_neig = mod_neig.get_permittivity();

          const libMesh::RealTensor& eps_average = (eps_elem + eps_neig) / 2.0;

          my_det_average = eps_average.det();

          a11 = eps_average(0,0), a12 = eps_average(0,1), a13 = eps_average(0,2),
          a21 = eps_average(1,0), a22 = eps_average(1,1), a23 = eps_average(1,2),
          a31 = eps_average(2,0), a32 = eps_average(2,1), a33 = eps_average(2,2);

          const libMesh::RealTensor inv_eps_average((a33*a22-a32*a23)/my_det_average, -(a33*a12-a32*a13)/my_det_average, (a23*a12-a22*a13)/my_det_average,
                                                 -(a33*a21-a31*a23)/my_det_average, (a33*a11-a31*a13)/my_det_average, -(a23*a11-a21*a13)/my_det_average,
                                                  (a32*a21-a31*a22)/my_det_average, -(a32*a11-a31*a12)/my_det_average, (a22*a11-a21*a12)/my_det_average);

          libMesh::VectorValue<Real> tmp = (eps_average * normal_face[qp]);

          for (unsigned int i = 0; i < n_dofs; i++)
          {
            lifting_d -= penalty / h_elem * (solution(dof_indices[i]) * phi_face[i][qp] - solution(neig_dof_indices[i]) * phi_face_neig[i][qp]) * tmp;
          }
          lifting_e +=inv_eps_average * lifting_d;
        }
      }
    }
  }

  if (do_displacement)
  {
    //TODO review the writing solution
    values[Displacement][0] = displac(0) / np + lifting_d(0);
    values[Displacement][1] = displac(1) / np + lifting_d(1);
    values[Displacement][2] = displac(2) / np + lifting_d(2);
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
    values[ElField][0] = field_cons(0) / np + lifting_e(0);
    values[ElField][1] = field_cons(1) / np + lifting_e(1);
    values[ElField][2] = field_cons(2) / np + lifting_e(2);
  }
}



void
Poisson::do_assemble(libMesh::EquationSystems&, const std::string&)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  ID dim = get_mesh().mesh_dimension();

  double mesh_units = get_mesh_units();

  libMesh::DofMap& dof_map =  system.get_dof_map(); // constant reference to this system's _dof_map.
  // _dof_map ---> data structure describing the relationship between nodes, variables, etc... and degrees of freedom.

  const unsigned int u_var = system.variable_number("u"); // he variable number associated with the user-specified variable named var.
  libMesh::FEType fe_type = dof_map.variable_type(u_var); //The finite element type for variable var.

  libMesh::AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  libMesh::AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type));
  libMesh::AutoPtr<FEBase> fe_face_neig(build_finite_element(dim, fe_type));

  //volume finite element
  libMesh::QGauss qrule(dim, FIRST);
  fe->attach_quadrature_rule(&qrule);

  // surface finite element
  libMesh::QGauss qface(dim - 1, CONSTANT);
  fe_face->attach_quadrature_rule(&qface);
  fe_face_neig->attach_quadrature_rule(&qface);

  const vector<Real>& JxW = fe->get_JxW();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& q_point =fe->get_xyz();

  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe_face->get_dphi();
  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& normal_face = fe_face->get_normals();
  const vector<Point>& qface_point = fe_face->get_xyz();

  const vector<vector<Real> >&  phi_face_neig = fe_face_neig->get_phi();
  const vector<vector<RealGradient> >& dphi_face_neig = fe_face_neig->get_dphi();

  libMesh::DenseMatrix<Number> Ke;
  libMesh::DenseVector<Number> Fe;
  libMesh::DenseMatrix<Number> Ken;

  vector<unsigned int> dof_indices;

  unsigned int penalty = 4;

  //iterate through the elements
  MeshBase::const_element_iterator el = get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el= get_mesh().active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    unsigned int n_sides = elem->n_sides();

    fe->reinit(elem);

    Ke.resize(n_dofs,n_dofs);
    Fe.resize(n_dofs);

    PoissonModel& mod = *get_bulk_model<PoissonModel>(elem);
    mod.set_element(elem);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      mod.set_point(q_point[qp]);

      mod.calculate(elem, q_point[qp]);

      const RealTensor& eps_elem = mod.get_permittivity();
      const RealVectorValue& pol = mod.get_polarization();
      double rho = Constants::e * mod.get_charge_density();

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
          Ke(i,j) += JxW[qp] * (eps_elem * dphi[i][qp] * dphi[j][qp]);

        Fe(i) += JxW[qp] * (rho * phi[i][qp] + pol * dphi[i][qp]) / Constants::e0;
      }
    }

    //the sides share
    for (unsigned int s = 0; s < n_sides; s++)
    {
      const Elem* elem_neig = elem->neighbor(s);

      fe_face->reinit(elem, s);

      PoissonBoundaryModel* mod_int = get_interface_model<PoissonBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod_int->calculate(elem, s, qface_point[qp]);

          double a, b, c;
          mod_int->get_coefficients(a, b, c);

          // we use a penalty approach here for its simplicity
          if ((b < 1e-12) && (b >= 0)) b = 1e-12;
          else if ((b > -1e-12) && (b<= 0)) b = -1e-12;

          a /= b;
          c /= b;

          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
              Ke(i,j) += a * JxW_face[qp] * (phi_face[i][qp] * phi_face[j][qp]);

            Fe(i) += c * JxW_face[qp] * phi_face[i][qp];
          }
        }
      }

      else if (elem_neig != NULL)
      {
        PoissonModel& mod_neig = *get_bulk_model<PoissonModel>(elem_neig);

        AutoPtr<Elem> elem_side(elem->build_side(s));
        const double h_elem = (elem->volume()/elem_side->volume()) * mesh_units;

        vector<Point> qface_point_neig;

        FEInterface::inverse_map(elem->dim(), fe_type, elem_neig, qface_point, qface_point_neig);

        std::vector<unsigned int> neig_dof_indices;
        dof_map.dof_indices (elem_neig,neig_dof_indices);
        const unsigned int n_dofs_neig = neig_dof_indices.size();

        fe_face_neig->reinit(elem_neig, &qface_point_neig);

        Ken.resize(n_dofs,n_dofs_neig);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod.calculate(elem, qface_point[qp]);
          const libMesh::RealTensor& eps_elem = mod.get_permittivity();

          mod_neig.calculate(elem_neig, qface_point[qp]);
          const libMesh::RealTensor& eps_neig = mod_neig.get_permittivity();

          const libMesh::RealTensor& eps_average = (eps_elem + eps_neig) / 2.0;

          libMesh::VectorValue<Real> tmp1 = (eps_elem * normal_face[qp]);
          libMesh::VectorValue<Real> tmp2 = (eps_average * normal_face[qp]);

          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
            {
              Ke(i,j) -= 0.5 * JxW_face[qp] * (phi_face[j][qp] * (tmp1 * dphi_face[i][qp]) +  phi_face[i][qp] * (tmp2 * dphi_face[j][qp]));//consistency
              Ke(i,j) += JxW_face[qp] * penalty / h_elem * phi_face[j][qp] * phi_face[i][qp] * tmp2 * normal_face[qp];  //stability
           }

            for (unsigned int j = 0; j < n_dofs_neig; j++)
            {
              Ken(i,j) += 0.5 * JxW_face[qp] * (phi_face_neig[j][qp] * (tmp1 * dphi_face[i][qp]) -  phi_face[i][qp] * (tmp2 * dphi_face_neig[j][qp]));//consistency
              Ken(i,j) -= JxW_face[qp] * penalty / h_elem * phi_face[i][qp] * phi_face_neig[j][qp] * tmp2 * normal_face[qp];  //stability
            }
          }
        }
        system.matrix->add_matrix(Ken, dof_indices, neig_dof_indices);
      }
    }
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }

  system.matrix->close();
  system.matrix->print_matlab("Ke.m");
  system.rhs->close();
  system.rhs->print_matlab("Fe.m");
}
