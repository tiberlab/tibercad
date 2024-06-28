// $Id: PVModule.C 5514 2024-05-03 15:22:47Z maufder $

#include "PVModule.h"
//#include "PVModuleModel.h"
//#include "PVModuleBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/vector_value.h"

// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace std;
using namespace libMesh;



PVModule::PVModule(const ModelOptions& options) :
  SimulationInterface(options),
  _my_assembly(this)
{
  // there's nothing to be done
}


PVModule::~PVModule(void)
{
  // there's nothing to be done
}


PVModule*
PVModule::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new PVModule(options);
}



void
PVModule::do_init(void)
{
  parse_options();

  // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  system.add_variable("Vtop", FIRST);
  system.add_variable("Vbot", FIRST);
  system.add_vector("currdens");
  system.attach_assemble_object(_my_assembly);
  system.init();
}


void
PVModule::parse_options(void)
{

  _spice = get_option("spice_executable", _spice);
 
}


void
PVModule::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(TopPotential, REAL, NODES, "V");
  declare_solution(BottomPotential, REAL, NODES, "V");
  declare_solution(CurrentDensity, REAL, NODES, "A/cm^2");
}


void
PVModule::do_solve(void)
{

  TiberLinearSystem& sys = get_equation_system<TiberLinearSystem>();


  //system.set_options(get_solver_options());
  //system.solve();

  sys.assemble();

  // write Spice netlist
  string netlist = get_output_directory() + "/" + get_name() + "_spice.net";
  ofstream of(netlist);
  of.close();

  // call ngspice
  Messages::info("calling Spice: " + _spice);
  string log = get_output_directory() + "/" + get_name() + "_spice.log";
  string outfile = get_output_directory() + "/" + get_name() + "_spice.dat";
  string cmdline = _spice + " -b -o " + log + " -r " + outfile + " " + netlist ;
  int ret = std::system(cmdline.c_str());

  if (ret == -1)
    throw(SolveFailedException("Could not run Spice."));

  // parse output and populate solution vectors
}


void
PVModule::do_print_info(void)
{
  Messages::info("PVModule: lumped element model for photovoltaic module"
                 " simulation.");
}


PhysicalModel*
PVModule::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  //return PVModuleModel::create(mat, options);
  return(nullptr);
}



PhysicalModel*
PVModule::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  //return PVModuleBoundaryModel::create(boundary, options);
  return(nullptr);
}



void
PVModule::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int v_top = system.variable_number("Vtop");
  const unsigned int v_bot = system.variable_number("Vbot");

  FEType fe_type = system.variable_type(v_top);
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dofs_top;
  vector<unsigned int> dofs_bot;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  //const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  //const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dofs_top, v_top);
  dof_map.dof_indices(elem, dofs_bot, v_bot);
  const unsigned int n_dofs = dofs_top.size();

  //PVModuleModel& mod = *get_bulk_model<PVModuleModel>(elem);


  for (unsigned int n = 0; n < np; n++)
  {
    double ut  = 0.0;
    double ub  = 0.0;

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      ut += phi[i][n] * solution(dofs_top[i]);
      ub += phi[i][n] * solution(dofs_bot[i]);
    }

    if (values.count(TopPotential))
      values[TopPotential][n] = ut;

    if (values.count(BottomPotential))
      values[BottomPotential][n] = ub;

  }
}



void
PVModule::assemble(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
 
  
  DofMap& dof_map =  system.get_dof_map();

  const unsigned int vtop = system.variable_number("Vtop");
  const unsigned int vbot = system.variable_number("Vbot");

  FEType fe_type = dof_map.variable_type(vtop);

  // the finite element
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type, true));
  unique_ptr<QBase> qrule(QBase::build(QTRAP, dim));
  fe->attach_quadrature_rule(qrule.get());
/*
  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  // the surface finite element
  unique_ptr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  unique_ptr<QBase> qface(QBase::build(myopts.quadrature_type, dim-1, myopts.integration_order));
  fe_face->attach_quadrature_rule(qface.get());

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

    PVModuleModel& mod = *get_bulk_model<PVModuleModel>(elem);

    mod.set_element(elem);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      mod.set_point(q_point[qp]);

      mod.calculate();

      const RealTensor& eps = mod.get_permittivity();
      // units of polarization ??????
      const RealVectorValue& pol = mod.get_polarization();
      double rho =  mod.get_charge_density() * Lambda;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
          Ke(i, j) += JxW[qp] * (dphi[i][qp] * (eps * dphi[j][qp]));

        Fe(i) += JxW[qp] * (rho * phi[i][qp] + pol * dphi[i][qp]);
      }

    }


    // the sides
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      PVModuleBoundaryModel* mod_int =
        get_interface_model<PVModuleBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        fe_face->reinit(elem, s);

        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
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
  system.matrix->close();
  //system.matrix->print_matlab("K.m");
*/
}
