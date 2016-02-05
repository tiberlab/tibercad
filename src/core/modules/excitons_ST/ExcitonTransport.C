
// $Id: ExcitonTransport.C 4192 2015-12-10 11:11:18Z maufder $

// module includes
#include "ExcitonTransport.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Scaling.h"
#include "Constants.h"
#include "ExcitonProperties.h"
#include "FiniteElement.h"
#include "Messages.h"

// libmesh includes
#include "node.h"
#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "equation_systems.h"
#include "TiberNonlinearSystem.h"
//#include "mesh_refinement.h"
//#include "error_vector.h"
//#include "kelly_error_estimator.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
// C++ includes
#include "TiberModule.h"

using namespace std;
using namespace ExcitonsDefs;



ExcitonTransport*
ExcitonTransport::_this;


ExcitonTransport::Options::Options(void)
  : mesh_refinement(false),
    max_refinement_steps(5),
    max_refinement_level(8),
    refine_fraction(0.7),
    coarsen_fraction(0.3),
    refinement_tolerance(1e-6),
    integration_order(libMeshEnums::FIFTH),
    coupling(FULL)
{
}

ExcitonTransport::Options::Options(const Options& rhs)
  : mesh_refinement(rhs.mesh_refinement),
    max_refinement_steps(rhs.max_refinement_steps),
    max_refinement_level(rhs.max_refinement_level),
    refine_fraction(rhs.refine_fraction),
    coarsen_fraction(rhs.coarsen_fraction),
    refinement_tolerance(rhs.refinement_tolerance),
    integration_order(rhs.integration_order),
    coupling(FULL)
{
}

ExcitonTransport::Options&
ExcitonTransport::Options::operator=(const Options& rhs)
{
  if (&rhs != this)
  {
    mesh_refinement = rhs.mesh_refinement;
    max_refinement_steps = rhs.max_refinement_steps;
    max_refinement_level = rhs.max_refinement_level;
    refine_fraction = rhs.refine_fraction;
    coarsen_fraction = rhs.coarsen_fraction;
    refinement_tolerance = rhs.refinement_tolerance;
    integration_order = rhs.integration_order;
    coupling = FULL;
  }
  return *this;
}




ExcitonTransport::ExcitonTransport(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true)
{
}


ExcitonTransport::~ExcitonTransport(void)
{
  cleanup_solver();
}




PhysicalModel*
ExcitonTransport::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  const string& modelname = options.get_option("model", "simple");

  ExcitonProperties* model =
    ExcitonProperties::create(modelname, mat, options);

  if (model == NULL)
    throw ModelErrorException(
        "ExcitonTransport: No such physical model: " + modelname);

  return model;
}


void
ExcitonTransport::compute_scaling(void)
{


  // we calculate in cm!
  double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  double x0 = -1;
  double phi0 = 1;
  double C0 = 1;

  const MeshBase& mesh = get_mesh();
  MeshBase::const_node_iterator it = mesh.nodes_begin();
  const MeshBase::const_node_iterator end = mesh.nodes_end();

  assert(it != end);

  const Node& n = **it;
  double xmin = n(0), ymin = n(1), zmin = n(2);
  double xmax = xmin, ymax = ymin, zmax = zmin;
  ++it;
  while (it != end)
  {
    const Node& n = **it;

    if (n(0) < xmin)
      xmin = n(0);
    else if (n(0) > xmax)
      xmax = n(0);

    if (n(1) < ymin)
      ymin = n(1);
    else if (n(1) > ymax)
      ymax = n(1);

    if (n(2) < zmin)
      zmin = n(2);
    else if (n(2) > zmax)
      zmax = n(2);

    ++it;
  }
  double x = xmax - xmin;
  double y = ymax - ymin;
  double z = zmax - zmin;

  x0 = (x > y) ? x : y;
  x0 = (x0 > z) ? x0 : z;

  phi0 = Constants::k_B * SimulationOptions::T;

  get_scaling().set_potential_scaling(phi0);
  get_scaling().set_length_scaling(x0 * mesh_units);
  get_scaling().set_density_scaling(C0);
}




void
ExcitonTransport::reset_solver(void)
{
  if (!_rebuild_eq_system)
  {
    //get_equation_systems().delete_system(get_equation_system_name());
    _rebuild_eq_system = true;
  }
}


void
ExcitonTransport::cleanup_solver(void)
{
  reset_solver();
}


void
ExcitonTransport::do_print_info(void)
{
  parse_options();

  Options& myopts = get_options();

  Messages::newline();

  ostringstream os;
  os << "solving for : ";
  if (myopts.coupling & SINGLET)
    os << "singlets";
  if (myopts.coupling & TRIPLET)
    os << " and triplets";

  os << endl;

  Messages::info(os.str());

}


void
ExcitonTransport::parse_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

  string coupling = opts.get_option("coupling", "");

  if (coupling == "singlet")
    myopts.coupling = SINGLET;
  else
    myopts.coupling = FULL;

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", (int) myopts.integration_order));

  myopts.mesh_refinement = opts.get_option("mesh_refinement", false);

  const unsigned int dim = get_mesh().mesh_dimension();
  if (dim == 1)
  {
    ModelOptions& solveropts = get_solver_options();
    string ksp_type(solveropts.get_option("ksp_type", "bcgsl"));
    if (ksp_type == "bcgsl")
      solveropts["ksp_type"] = "bcgs";
  }

}



void
ExcitonTransport::do_init(void)
{
  if (!_rebuild_eq_system) return;

  _device = &get_environment().get_device();

  ModelOptions::submodel_iterator linit(
      get_solver_options().submodels_begin("linear_solver"));

  if (linit == get_solver_options().submodels_end("linear_solver"))
  {
    get_solver_options().add_submodel("linear_solver", ModelOptions());
    linit = get_solver_options().submodels_begin("linear_solver");
  }

  ModelOptions& linopts = linit->second;

  // default is bcgsl
  if (!linopts.find_option("method"))
    linopts["method"] = "bcgsl";

  // in 1D bcgs seems to work better than bcgsl
  const unsigned int dim = get_mesh().mesh_dimension();
  if ((dim == 1) && (linopts["method"] == "bcgsl"))
    linopts["method"] = "bcgs";

  if (!linopts.find_option("preconditioner"))
    if (dim < 3)
      linopts["preconditioner"] = "lu";
    else
      linopts["preconditioner"] = "ilu";

  if (linopts.get_option("absolute_tolerance", -1.0) < 0)
    linopts["absolute_tolerance"] = "1e-15";


  ModelOptions& solveropts = get_solver_options();
  if (solveropts.get_option("absolute_tolerance", -1.0) < 0)
    solveropts["absolute_tolerance"] = "1e-15";
  EquationSystems& equation_systems = get_equation_systems();

  // create the exciton continuity equation
  create_equation_system("nonlinear");
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);


  system.add_variable("S", libMeshEnums::FIRST);
  system.add_variable("T", libMeshEnums::FIRST);

  // a local scaling
  //system.add_vector("local_scaling");

  system.attach_assembly_routine(assemble);

  //
  // initialize the newly created system
  system.init();

  // compute scaling factors
  compute_scaling();


  _rebuild_eq_system = false;

  //set_initial_guess(get_option("x_fermi_guess", 0.0));

}





void
ExcitonTransport::set_initial_guess(double guess)
{
  assert(_rebuild_eq_system == false);

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  system.get_solution_vector().zero();
  //system.get_solution_vector().add(guess);
}



void
ExcitonTransport::do_solve(void)
{
  assert(_rebuild_eq_system == false);

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;

  parse_options();


  EquationSystems& equation_systems = get_equation_systems();

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  NumericVector<Number>& solution = system.get_solution_vector();

  // set the solver parameters (they could have change since we made
  // the first calculation)
  system.set_options(get_solver_options());

  int coupling = get_options().coupling;

  //Utils::Timer tt;

  try
  {
    system.solve();
  }
  catch (SolverException& e)
  {
    string msg("ExcitonTransport: solve failed (" +
      string(e.what()) + ")");
    throw SolveFailedException(msg);
  }


  _n_nonlinear_iterations = system.n_nonlinear_iterations();
  _final_residual = system.final_residual_norm();

}

double
ExcitonTransport::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_density_scaling();
}


void
ExcitonTransport::assemble(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  switch (_this->_options.coupling)
  {
    case (SINGLET):
      _this->do_assembly<SINGLET>(x, residual, jacobian);
      break;
    default:
      _this->do_assembly<FULL>(x, residual, jacobian);
      break;
  }

}

template <int coupling>
void
ExcitonTransport::do_assembly(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  // references for nicer code
  const MeshBase& mesh = get_mesh();

  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  //NumericVector<Number>& locscal = system.get_vector("local_scaling");

  const unsigned int dim = mesh.mesh_dimension();

  const Device& device = *_device;
  const Options& params = get_options();
  Options& options = get_options();

  //
  // some scaling stuff...
  //
  // NOTE: the mesh and all paramters were not explicitly scaled, so
  //       we have to treat scaling by explicit division/multiplication
  //
  // the scaling parameters
  const Scaling& scaling = get_scaling();
  // the scaling parameter for the poisson eq.
  // The factor 1e-2 comes from the fact, that we are calculating in cm!
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double C0 = scaling.get_density_scaling();
  // scaling for recombination rates
  double R0 = C0 / scaling.get_time_scaling();



  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int es_var = system.variable_number("S");
  const unsigned int et_var = system.variable_number("T");

  FEType fe_type = system.variable_type(es_var);

  // the finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, params.integration_order);
  fe->attach_quadrature_rule(&qrule);

  // the finite element for boundary integration
  //AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  //libMeshEnums::Order integration_order;
  //if (dim == 1)
  //  integration_order = libMeshEnums::CONSTANT;
  //else
  //  integration_order = params.integration_order;
  //
  //QGauss qface(dim - 1, integration_order);
  //fe_face->attach_quadrature_rule(&qface);


  // references to cell-specific data that will be used to
  // assemble the system.
  // Data will be given for each quadrature point.
  //
  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe->get_JxW();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();
  //
  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  //
  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // the system matrix (will hold also element jacobian contribution)
  DenseMatrix<Number> Ke;
  // the system rhs (will hold also element rhs contribution)
  DenseVector<Number> Fe;
  // the local solution
  DenseVector<Number> X;

  DenseSubMatrix<Number>
    Kss(Ke), Kst(Ke),
    Kts(Ke), Ktt(Ke);

  DenseSubVector<Number>
    Fs(Fe),
    Ft(Fe);

  DenseSubVector<Number>
    Xs(X),
    Xt(X);

  // the DOF indices
  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_s;
  vector<unsigned int> dof_indices_t;


  // zero out residual and jacobian !! IMPORTANT !!
  if (residual != NULL)
    residual->zero();
  if (jacobian != NULL)
    jacobian->zero();

  MeshBase::const_element_iterator el =
                                  mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_local_elements_end();


  // loop over all active elements
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_s, es_var);
    dof_map.dof_indices(elem, dof_indices_t, et_var);

    unsigned int n_dofs     = dof_indices_s.size();
    unsigned int n_dofs_tot = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs_tot, n_dofs_tot);
    Fe.resize(n_dofs_tot);
    X.resize(n_dofs_tot);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);

    // Reposition the submatrices according to this scheme:
    //
    //         -           -          -  -
    //   Ke = | Kss Kst |   Fe = | Fs |
    //        | Kts Ktt |;       | Ft |
    //
    Kss.reposition(0, 0, n_dofs, n_dofs);
    Kst.reposition(0, n_dofs, n_dofs, n_dofs);
    //
    Kts.reposition(n_dofs, 0, n_dofs, n_dofs);
    Ktt.reposition(n_dofs, n_dofs, n_dofs, n_dofs);
    //
    Fs.reposition(0, n_dofs);
    Ft.reposition(n_dofs, n_dofs);
    //
    Xs.reposition(0, n_dofs);
    Xt.reposition(n_dofs, n_dofs);



    ExcitonProperties* em =
      dynamic_cast<ExcitonProperties*>(
          device.get_material(subdomain)->get_model(get_id()));


    assert(em != NULL);
    em->reinit(elem);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the exciton density at the quadrature point
      double s = 0.0;
      double t = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        s += phi[i][qp] * Xs(i);
        t += phi[i][qp] * Xt(i);
      }

      s *= C0;
      t *= C0;

      em->set_coordinates(q_point[qp]);
      em->set_density(s, t);
      em->calculate_diffusion();
      em->calculate_net_recombination_rate();

      double Rs = em->get_s_net_recombination_rate();
      double Rt = em->get_t_net_recombination_rate();
      double s_diffusion_coeff = em->get_s_diffusion();
      double t_diffusion_coeff = em->get_t_diffusion();


      //
      // The residual looks like this:
      //
      //      r_i = Ke_ij*X_j - Fe_i
      //
      // The jacobian looks like this:
      //
      //      J_ij = dr_i/dX_j
      //
      //           = Ke_ij + dKe_il/dX_j * X_l - dFe_i/dX_j
      //

      // the jacobian x weight x scaling
      double J = JxW[qp];

      //
      // First we will build the system matrix Ke_ij
      //

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          double laplace =
            J * (dphi[i][qp] * dphi[j][qp]);

            if (coupling & SINGLET)
              Kss(i,j) += s_diffusion_coeff * laplace;
            if (coupling & TRIPLET)
              Ktt(i,j) += t_diffusion_coeff * laplace;
        }

        if (!(coupling & SINGLET))
          Kss(i,i) += 1;
        if (!(coupling & TRIPLET))
          Ktt(i,i) += 1;
      }

      //
      // for jacobian compute the other contributions
      //
      if (jacobian != NULL)
      {
        double dRs = C0 / R0 * em->get_s_net_recombination_rate_derivative();
        double dRt = C0 / R0 * em->get_t_net_recombination_rate_derivative();
        double disc = C0 / R0 * em->get_isc_rate_derivative();

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (coupling & SINGLET)
              Kss(i,j) += dRs * phi_i_x_phi_j;
            if (coupling & TRIPLET)
              Ktt(i,j) += dRt * phi_i_x_phi_j;
              Kts(i,j) += disc * phi_i_x_phi_j;
          }
        }
      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {

        // net recombination rate
        double J_x_Rs = J * Rs / R0 ;
        double J_x_Rt = J * Rt / R0 ;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double net_recomb_s = J_x_Rs * phi[i][qp];
          double net_recomb_t = J_x_Rt * phi[i][qp];

          if (coupling & SINGLET)
            Fs(i) += net_recomb_s;
          if (coupling & TRIPLET)
            Ft(i) += net_recomb_t;

        }
      }

    }//end loop over qp

    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);


    if (residual != NULL)
    {
      for (unsigned int i = 0; i < n_dofs_tot; i++)
        for (unsigned int j = 0; j < n_dofs_tot; j++)
          Fe(i) += Ke(i,j) * x(dof_indices[j]);

      residual->add_vector(Fe, dof_indices);
    }
    else
      jacobian->add_matrix(Ke, dof_indices);

  } // end loop over elements

  if (residual != NULL)
  {
    residual->close();
    //residual->print_matlab("res.m");
  }

  if (jacobian != NULL)
  {
    jacobian->close();
    //jacobian->print_matlab("jac.m");
  }

}



void ExcitonTransport::get_solution_secure(const Elem* elem, std::map<ID, std::vector<double> >& solutions, const std::vector<Point>& points) {
    unsigned int np = points.size();

    TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

    const NumericVector<Number>& ddsol = system->get_solution_vector();

    const unsigned int dim = get_mesh().mesh_dimension();

    const DofMap& dof_map = system->get_dof_map();

    const unsigned int s_var = system->variable_number("S");
    const unsigned int t_var = system->variable_number("T");

    FEType fe_type = system->variable_type(s_var);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

    vector<unsigned int> dof_indices_s;
    vector<unsigned int> dof_indices_t;

    // element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();

    // element shape function gradients
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();

    fe->reinit(elem, &points);

    const vector<Point>& q_points = fe->get_xyz();

    dof_map.dof_indices(elem, dof_indices_s, s_var);
    dof_map.dof_indices(elem, dof_indices_t, t_var);

    const unsigned int n_dofs = dof_indices_s.size();


    ID subdomain = elem->subdomain_id();

    ExcitonProperties* excitonmodel =
      dynamic_cast<ExcitonProperties*>(
          _device->get_material(subdomain)->get_model(get_id()));

    assert(excitonmodel != NULL);

    excitonmodel->reinit(elem);

    // the scaling parameters to scale back the result
    double C0 = get_scaling().get_density_scaling();

    for (unsigned int n = 0; n < np; n++) {

      // do interpolation
      double s = 0;
      double t = 0;
      Point grad_s = 0.0;
      Point grad_t = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++) {
        s  += phi[i][n] * ddsol(dof_indices_s[i]);
        t  += phi[i][n] * ddsol(dof_indices_t[i]);
        grad_s  += dphi[i][n] * ddsol(dof_indices_s[i]);
        grad_t  += dphi[i][n] * ddsol(dof_indices_t[i]);
      }

      // scale the potential back
      s  *= C0;
      t  *= C0;
      grad_s *= C0;
      grad_t *= C0;

      // prepare for calculating local properties
      excitonmodel->set_coordinates(q_points[n]);

      double T_lat = excitonmodel->get_lattice_temperature();
      // all are at lattice temperature
      excitonmodel->set_carrier_temperature(T_lat);

      excitonmodel->set_density(s, t);

      double Ds = excitonmodel->get_s_diffusion();
      double Dt = excitonmodel->get_t_diffusion();


      if (solutions.count(SDENSITY)) {
        solutions[SDENSITY][n] = s;
      }

      if (solutions.count(SDIFFUSION)) {
        solutions[SDIFFUSION][n] = excitonmodel->get_s_diffusion();
      }

      if (solutions.count(SJ)) {
        solutions[SJ][3*n] = -Ds * grad_s(0);
        solutions[SJ][3*n + 1] = -Ds * grad_s(1);
        solutions[SJ][3*n + 2] = -Ds * grad_s(2);
      }

      if (solutions.count(SRDISS)) {
        solutions[SRDISS][n] = excitonmodel->get_s_dissociation_rate();
      }

      if (solutions.count(SRRAD)) {
        solutions[SRRAD][n] = excitonmodel->get_s_radiative_recombination_rate();
      }

      if (solutions.count(SRNONRAD)) {
        solutions[SRNONRAD][n] = excitonmodel->get_s_nonradiative_recombination_rate();
      }

      if (solutions.count(SGEN)) {
        solutions[SGEN][n] = excitonmodel->get_s_generation_rate();
      }

      if (solutions.count(SNETRECOMB)) {
        excitonmodel->calculate_net_recombination_rate();
        solutions[SNETRECOMB][n] = excitonmodel->get_s_net_recombination_rate();
      }

      if (solutions.count(SRADPOWER)) {
        solutions[SRADPOWER][n] =
        excitonmodel->get_exciton_energy() * excitonmodel->get_s_density() /
        excitonmodel->get_s_radiative_recombination_rate();
      }



      if (solutions.count(TDENSITY)) {
        solutions[TDENSITY][n] = t;
      }

      if (solutions.count(TDIFFUSION)) {
        solutions[TDIFFUSION][n] = excitonmodel->get_t_diffusion();
      }

      if (solutions.count(TJ)) {
        solutions[TJ][3*n] = -Dt * grad_t(0);
        solutions[TJ][3*n + 1] = -Dt * grad_t(1);
        solutions[TJ][3*n + 2] = -Dt * grad_t(2);
      }

      if (solutions.count(TRDISS)) {
        solutions[TRDISS][n] = excitonmodel->get_t_dissociation_rate();
      }

      if (solutions.count(TRRAD)) {
        solutions[TRRAD][n] = excitonmodel->get_t_radiative_recombination_rate();
      }

      if (solutions.count(TRNONRAD)) {
        solutions[TRNONRAD][n] = excitonmodel->get_t_nonradiative_recombination_rate();
      }

      if (solutions.count(TGEN)) {
        solutions[TGEN][n] = excitonmodel->get_t_generation_rate();
      }

      if (solutions.count(TNETRECOMB)) {
        excitonmodel->calculate_net_recombination_rate();
        solutions[TNETRECOMB][n] = excitonmodel->get_t_net_recombination_rate();
      }

      if (solutions.count(TRADPOWER)) {
        solutions[TRADPOWER][n] =
        excitonmodel->get_exciton_energy() * excitonmodel->get_t_density() /
        excitonmodel->get_t_radiative_recombination_rate();
      }

     if (solutions.count(RDISS)) {
        solutions[RDISS][n] = excitonmodel->get_s_dissociation_rate() + excitonmodel->get_t_dissociation_rate();
      }

    }

}



void
ExcitonTransport::do_setup_solution_variables(void) {

  declare_solution_ext("s_dens", SDENSITY, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_current", SJ, SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_dissociation", SRDISS, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_rad_recombination", SRRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_nonrad_recombination", SRNONRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_generation", SGEN, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_net_recombination", SNETRECOMB, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("s_rad_power", SRADPOWER, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");

  declare_solution_ext("t_dens", TDENSITY, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_current", TJ, SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_dissociation", TRDISS, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_rad_recombination", TRRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_nonrad_recombination", TRNONRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_generation", TGEN, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_net_recombination", TNETRECOMB, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("t_rad_power", TRADPOWER, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");

  declare_solution_ext("dissociation", RDISS, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
}

void ExcitonTransport::get_solution_secure(std::map<ID, std::vector<double> >& solutions) {

}

