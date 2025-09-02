// $Id$

#include "WaterIngress.h"
#include "WIModel.h"
#include "WIBoundaryModel.h"
#include "WIUtils.h"
#include "SimulationOptions.h"
#include "TiberTransientSystem.h"
#include "Messages.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;
using namespace libMesh;



WaterIngress::WaterIngress(const ModelOptions& options) :
  SimulationInterface(options),
  _my_assembly(this)
{
  // there's nothing to be done
}


WaterIngress::~WaterIngress(void)
{
  // there's nothing to be done
}


WaterIngress*
WaterIngress::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new WaterIngress(options);
}



void
WaterIngress::do_init(void)
{
  parse_options();

  // create a linear equation system 
  create_equation_system("transient");

  // get the reference to it
  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  // we use the partial pressure as variable, as this is the quantity
  // that is continuous
  // The concentration then is given by Henry's law as c = S*p
  system.add_variable("p", FIRST, &(this->get_region_ids()));
  system.attach_assemble_object(_my_assembly);
  system.init();
  //system.get_solution_vector().zero();

}


void
WaterIngress::parse_options(void)
{
  // no options right now 
  _cell_temp = get_option("cell_temperature", SimulationOptions::temperature);
}


void
WaterIngress::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(PartialPressure, REAL, NODES, "Pa");
  declare_solution(Concentration, REAL, NODES, "g/cm^3");
  declare_solution(RelativeHumidity, REAL, NODES, "%");
  declare_solution(Flux, VECTOR, CELL, "g/cm^2/s");
  declare_solution(Solubility, REAL, NODES, "g/cm^3/Pa");
  declare_solution(Diffusivity, REAL, NODES, "cm^2/s");
}


void
WaterIngress::do_solve(void)
{
  double current_time = TiberCad::get_global_time();

  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  system.set_options(get_solver_options());
  system.set_target_time(current_time);
  system.solve();

}


void
WaterIngress::do_print_info(void)
{
  Messages::info("Simulation of water ingress based on Fick's laws");
}


PhysicalModel*
WaterIngress::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  return WIModel::create(mat, options);
}



PhysicalModel*
WaterIngress::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return WIBoundaryModel::create(boundary, options);
}



void
WaterIngress::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("p");

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
  RealGradient flux(0);
  WIModel& mod = *get_bulk_model<WIModel>(elem);

  mod.calculate(elem, elem->vertex_average());
  
  // parameters are given in units of m
  double S = 1e-6 * mod.get_solubility();
  double D = 1e4 * mod.get_diffusivity();

 
  for (unsigned int n = 0; n < np; n++)
  {
    double p  = 0.0;
    RealGradient grad(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      p += phi[i][n] * solution(dof_indices[i]);
      grad -= dphi[i][n] * solution(dof_indices[i]);
    }

    flux += S*D*grad;

    // it might happen that numerical noise needs
    // to small negative values. We assure p >= 0
    if (p < 0) p = 0;

    if (values.count(PartialPressure))
      values[PartialPressure][n] = p;

    if (values.count(Concentration))
      values[Concentration][n] = S*p;

    if (values.count(RelativeHumidity))
    {
      double psat = WIUtils::goff_gratch(_cell_temp);
      values[RelativeHumidity][n] = 100 * p / psat;
    }

    if (values.count(Solubility))
      values[Solubility][n] = S;

    if (values.count(Diffusivity))
      values[Diffusivity][n] = D;

  }


  // this is a very primitive way of estimating the mean value.
  if (values.count(Flux))
  {
    values[Flux][0] = flux(0) / np;
    values[Flux][1] = flux(1) / np;
    values[Flux][2] = flux(2) / np;
  }

}



void
WaterIngress::assemble(void)
{
  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  // this vector contains the diagonal of the (lumped) A matrix
  NumericVector<libMesh::Number>& t_weight = system.get_vector("t_weight");
  t_weight.zero();

  // the time relevant for assembly
  double current_time = system.time;

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  // NB: Tibercad by default uses length-scale in meters 
  //     This means that FEM derivatives d/dx are in 1/m
  //     To change this behavior it is necessary to define a different 'scaling'
  //     For instance if we want to use mesh_units in the assembly we need to:
  //     1. set the scaling to mesh units:
  //        get_scaling().set_length_scaling(get_mesh_units());
  //     2. use  build_finite_element(dim, fe_type, true)  
  //                                                ^ false is the default  
  //get_scaling().set_length_scaling(100 * get_mesh_units());
  get_scaling().set_calc_mesh_units(100 * get_mesh_units());

  // a set with all Dirichlet DoF indices
  std::set<dof_id_type> dirichlet_dofs;

  double penalty_value = 1e-12;

  DofMap& dof_map =  system.get_dof_map();

  const unsigned int uvar = system.variable_number("p");

  FEType fe_type = dof_map.variable_type(uvar);

  // the finite element
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type, true));
  unique_ptr<QBase> qrule(QBase::build(QTRAP, dim, FIRST));
  fe->attach_quadrature_rule(qrule.get());

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  // the surface finite element
  unique_ptr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  unique_ptr<QBase> qface(QBase::build(QTRAP, dim-1, FIRST));
  fe_face->attach_quadrature_rule(qface.get());

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();

  vector<unsigned int> dof_indices;

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  DenseVector<Number> A;
  DenseVector<Number> sol;

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
    A.resize(n_dofs);
    sol.resize(n_dofs);

    // get the current solution
    dof_map.extract_local_vector(solution, dof_indices, sol);

    WIModel& mod = *get_bulk_model<WIModel>(elem);

    mod.calculate(elem, elem->vertex_average());

    double S = 1e-6 * mod.get_solubility();
    double D = 1e4 * mod.get_diffusivity();

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          Ke(i, j) += JxW[qp] * S * D * (dphi[i][qp] * dphi[j][qp]);

          // we use mass lumping here
          A(i) += JxW[qp] * S * (phi[i][qp] * phi[j][qp]);
        }
      }

    }

    // the sides
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      WIBoundaryModel* mod_int =
        get_interface_model<WIBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        fe_face->reinit(elem, s);

        bool is_dirichlet = false;

        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {
          mod_int->calculate(elem, s, qface_point[qp]);

          double a, b, c;
          mod_int->get_coefficients(a, b, c);

          // we use a penalty approach here for its simplicity
          if ((b < penalty_value) && (b >= 0))
          {
            b = penalty_value;
            is_dirichlet = true;
          }
          else if ((b > -penalty_value) && (b<= 0))
          {
            b = -penalty_value;
            is_dirichlet = true;
          }

          a /= b;
          c /= b;

          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
              Ke(i, j) += a * JxW_face[qp] * (phi_face[i][qp] * phi_face[j][qp]);

            Fe(i) += c * JxW_face[qp] * phi_face[i][qp];
            
          }
        }

        if (is_dirichlet)
        {
          for (unsigned int i = 0; i < elem->n_nodes(); ++i)
          {
            if (elem->is_node_on_side(i, s))
            {
              dirichlet_dofs.insert(dof_indices[i]);
            }
          }
        }
      }
    }


    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
    t_weight.add_vector(A, dof_indices);

  }
  system.matrix->close();
  system.rhs->close();
  t_weight.close();
  //system.matrix->print_matlab("K_wi.m");
  //system.rhs->print_matlab("F_wi.m");
  //t_weight.print_matlab("A_wi.m");

  // let the system know which are Dirichlet DoFs
  system.set_dirichlet_dofs(dirichlet_dofs);

  // for Dirichlet DoFs we set A to zero, so the corresponding
  // equations become algebraic. This zero might be used by the
  // transient solver, and it is more elegant for fulfilling
  // essential boundy conditions
  for (auto& i : dirichlet_dofs)
  {
    //double value = system.rhs->el(i) / (*system.matrix)(i, i);
    //system.solution->set(i, value);
    t_weight.set(i, 0.0);
  }

  t_weight.close();
}
