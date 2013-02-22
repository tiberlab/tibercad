
// $Id$

// module includes
#include "ExcitonTransport.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Scaling.h"
#include "Constants.h"
#include "ExcitonProperties.h"
#include "FiniteElement.h"

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
#include "TiberNonlinearScalarSolver.h"
// C++ includes

using namespace std;


//
// Module interface
//

TIBER_MODULE(ExcitonTransport, MODULE_NAME)



ExcitonTransport*
ExcitonTransport::_this;


ExcitonTransport::Options::Options(void)
  : mesh_refinement(false),
    max_refinement_steps(5),
    max_refinement_level(8),
    refine_fraction(0.7),
    coarsen_fraction(0.3),
    refinement_tolerance(1e-6),
    integration_order(libMeshEnums::FIFTH)
{
}

ExcitonTransport::Options::Options(const Options& rhs)
  : mesh_refinement(rhs.mesh_refinement),
    max_refinement_steps(rhs.max_refinement_steps),
    max_refinement_level(rhs.max_refinement_level),
    refine_fraction(rhs.refine_fraction),
    coarsen_fraction(rhs.coarsen_fraction),
    refinement_tolerance(rhs.refinement_tolerance),
    integration_order(rhs.integration_order)
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
  double mu0 = 1;
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
  get_scaling().set_mobility_scaling(mu0);
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
ExcitonTransport::parse_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

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

  if (polaritons) {
    // TODO :
    linopts["method"] = "bcgs"; // test
  }

  ModelOptions& solveropts = get_solver_options();
  if (solveropts.get_option("absolute_tolerance", -1.0) < 0)
    solveropts["absolute_tolerance"] = "1e-15";
  EquationSystems& equation_systems = get_equation_systems();

  // create the exciton continuity equation
  create_equation_system("nonlinear");
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);


  system.add_variable("fermi_x", libMeshEnums::FIRST);

  // a local scaling
  system.add_vector("local_scaling");

  system.attach_assembly_routine(assemble);

  //
  // initialize the newly created system
  system.init();

  // compute scaling factors
  compute_scaling();


  _rebuild_eq_system = false;

  set_initial_guess(get_option("x_fermi_guess", 0.0));

  if (polaritons) {
    scalarSolver.init(1.0);
    scalarSolver.attach_assembly_routine(assemble_DRSystem);
  }
}





void
ExcitonTransport::set_initial_guess(double guess)
{
  assert(_rebuild_eq_system == false);

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  system.get_solution_vector().zero();
  system.get_solution_vector().add(guess);
}



void
ExcitonTransport::build_local_scaling(void)
{
  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();
  
  
  const NumericVector<Number>& solution = *(system->solution);
  NumericVector<Number>& locscal = system->get_vector("local_scaling");

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const Options& params = get_options();

  // the scaling parameters to scale back the result
  const Scaling& scaling = get_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double C0 = scaling.get_density_scaling();

  const unsigned int var = system->variable_number("fermi_x");
  
  vector<unsigned int> dof_indices;

  FEType fe_type = system->variable_type(var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, params.integration_order);
  fe->attach_quadrature_rule(&qrule);


  // Jacobian * quadrature weight at each integration point.   
  const vector<Real>& JxW = fe->get_JxW();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();
  //
  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  //
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  locscal.zero();


  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices, var);

    unsigned int n_dofs     = dof_indices.size();

    ExcitonProperties* sc =
      dynamic_cast<ExcitonProperties*>(
          device.get_material(subdomain)->get_model(get_id()));
    assert(sc != NULL); 

    sc->reinit(elem);

    fe->reinit(elem);


    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real ex  = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++)
        ex  += phi[i][qp] * solution(dof_indices[i]);



      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      double T_lat = sc->get_lattice_temperature();
      // all are at lattice temperature
      sc->set_carrier_temperature(T_lat);

      sc->set_effective_potential(ex);
      sc->calculate_density();
      sc->calculate_net_recombination_rate();


      double sigma_x = JxW[qp] * (sc->get_mobility() * sc->get_density() + 0e-12);

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //if (sigma_x > 1)
        //  locscal.add(dof_indices[i], sigma_x * (dphi[i][qp] * dphi[i][qp]));
        //else
          locscal.set(dof_indices[i], 1.0);
      }


    } // end loop over quadrature points
  } // end loop over elements

  locscal.close();
}



void
ExcitonTransport::do_solve(void)
{
  assert(_rebuild_eq_system == false);

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;

  parse_options();

  build_local_scaling();

  EquationSystems& equation_systems = get_equation_systems();

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  NumericVector<Number>& solution = system.get_solution_vector();

  // set the solver parameters (they could have change since we made
  // the first calculation)
  system.set_options(get_solver_options());

  //Utils::Timer tt;
  calc_excpol_integral(solution);
  //std::cout << "Excpol int time: "  << tt.elapsed_string() << "\n";

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

  calc_excpol_integral(solution);
}

double
ExcitonTransport::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_potential_scaling();
}

void
ExcitonTransport::do_assembly(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  // references for nicer code
  const MeshBase& mesh = get_mesh();

  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  NumericVector<Number>& locscal = system.get_vector("local_scaling");

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
  const double mu0 = scaling.get_mobility_scaling();
  // scaling for recombination rates
  double R0 = C0 / scaling.get_time_scaling();


  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int ex_var = system.variable_number("fermi_x");

  FEType fe_type = system.variable_type(ex_var);

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

  // the DOF indices
  vector<unsigned int> dof_indices;


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
    unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);
    X.resize(n_dofs);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);

    ExcitonProperties* em =
      dynamic_cast<ExcitonProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(em != NULL);

    em->reinit(elem);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the exciton electro-chemical potential at the quadrature point
      double ex = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++)
        ex += phi[i][qp] * X(i);

      em->set_coordinates(q_point[qp]);
      em->set_effective_potential(phi0 * ex);
      em->calculate_density();
      em->calculate_mobility();
      em->calculate_net_recombination_rate();

      double x = em->get_density();
      double Rx = em->get_net_recombination_rate();
      double mux = em->get_mobility();


      // NOTE: sigma_x = mu_x * x is the exciton conductivity
      double sigma_x = (mux * x + 0e-12) / (mu0 * C0);


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

            Ke(i,j) += sigma_x * laplace;
        }
      }

      //
      // for jacobian compute the other contributions
      //
      if (jacobian != NULL)
      {
        double dRx = phi0 / R0 * em->get_net_recombination_rate_derivative();
        if (Rx == 0.0)
          dRx = 0.0;

        double dsigma_x = J * phi0 / C0 * mux * em->get_density_derivative() ;
          

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part

            double dsigma_x_phi = dsigma_x * phi[j][qp] ;

            for (unsigned int k = 0; k < n_dofs; k++)
            {
              double laplace = (dphi[i][qp] * dphi[k][qp]);

              double elem_contrib =
                dsigma_x_phi * laplace * X(k);

              Ke(i,j) += elem_contrib / locscal(dof_indices[i]);

            }

            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            Ke(i,j) += dRx * phi_i_x_phi_j / locscal(dof_indices[i]);
          }
        }
      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {

        // net recombination rate
        double J_x_Rx = J * Rx / R0 ;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double net_recomb_x = J_x_Rx * phi[i][qp];

          Fe(i) += net_recomb_x / locscal(dof_indices[i]);
        }
      }

    }//end loop over qp

    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);


    if (residual != NULL)
    {
      for (unsigned int i = 0; i < n_dofs; i++)
        for (unsigned int j = 0; j < n_dofs; j++)
          Fe(i) += Ke(i,j) * x(dof_indices[j]);

      residual->add_vector(Fe, dof_indices);
    }
    else
      jacobian->add_matrix(Ke, dof_indices);

  } // end loop over elements

  if (residual != NULL)
    residual->close();

  if (jacobian != NULL)
    jacobian->close();

}



void ExcitonTransport::get_solution_secure(const Elem* elem, std::map<ID, std::vector<double> >& solutions, const std::vector<Point>& points) {
    unsigned int np = points.size();

    TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

    const NumericVector<Number>& ddsol = *(system->solution);

    const unsigned int dim = get_mesh().mesh_dimension();

    const DofMap& dof_map = system->get_dof_map();

    const unsigned int u_var = system->variable_number("fermi_x");

    FEType fe_type = system->variable_type(u_var);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

    vector<unsigned int> dof_indices_u;

    // element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();

    // element shape function gradients
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();

    fe->reinit(elem, &points);

    const vector<Point>& q_points = fe->get_xyz();

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    const unsigned int n_dofs = dof_indices_u.size();


    ID subdomain = elem->subdomain_id();

    ExcitonProperties* excitonmodel =
      dynamic_cast<ExcitonProperties*>(
          _device->get_material(subdomain)->get_model(get_id()));

    assert(excitonmodel != NULL);

    excitonmodel->reinit(elem);

    // the scaling parameters to scale back the result
    double phi0 = get_scaling().get_potential_scaling();

    for (unsigned int n = 0; n < np; n++) {

      // do interpolation
      double ex = 0;
      Point grad_x = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++) {
        ex  += phi[i][n] * ddsol(dof_indices_u[i]);
        grad_x  += dphi[i][n] * ddsol(dof_indices_u[i]);
      }

      // scale the potential back
      ex  *= phi0;
      grad_x *= phi0;

      // prepare for calculating local properties
      excitonmodel->set_coordinates(q_points[n]);

      double T_lat = excitonmodel->get_lattice_temperature();
      // all are at lattice temperature
      excitonmodel->set_carrier_temperature(T_lat);

      excitonmodel->set_effective_potential(ex);

      excitonmodel->calculate_density();

      if (solutions.count(RDISS)) {
        solutions[RDISS][n] = excitonmodel->get_dissociation_rate();
      }

      double sigma = excitonmodel->get_real_density() * excitonmodel->get_mobility();

      if (solutions.count(CHEMPOT)) {
        solutions[CHEMPOT][n] = ex;
      }

      if (solutions.count(XSIGMA)) {
        solutions[XSIGMA][n] = sigma;
      }

      if (solutions.count(XMOBILITY)) {
        solutions[XMOBILITY][n] = excitonmodel->get_mobility();
      }

      if (solutions.count(XDENSITY)) {
        solutions[XDENSITY][n] = excitonmodel->get_real_density();
	excpolprops.max_density = std::max(excpolprops.max_density, excitonmodel->get_real_density());
      }

      if (solutions.count(RDISS)) {
        solutions[RDISS][n] = excitonmodel->get_dissociation_rate();
      }

      if (solutions.count(RRAD)) {
        solutions[RRAD][n] = excitonmodel->get_radiative_recombination_rate();
      }

      if (solutions.count(RNONRAD)) {
        solutions[RNONRAD][n] = excitonmodel->get_nonradiative_recombination_rate();
      }

      if (solutions.count(EXCEXCPOLARITON)) {
        solutions[EXCEXCPOLARITON][n] = excitonmodel->get_exc_exc_scattering();
      }

      if (solutions.count(EXCPHONONPOLARITON)) {
        solutions[EXCPHONONPOLARITON][n] = excitonmodel->get_exc_photon_scattering();
      }

      if (solutions.count(NETRECOMB)) {
        excitonmodel->calculate_net_recombination_rate();
        solutions[NETRECOMB][n] = excitonmodel->get_real_net_recombination_rate();
      }

      if (solutions.count(XGEN)) {
        solutions[XGEN][n] = excitonmodel->get_generation_rate();
      }

      if (solutions.count(J)) {
        solutions[J][3*n] = -sigma * grad_x(0);
        solutions[J][3*n + 1] = -sigma * grad_x(1);
        solutions[J][3*n + 2] = -sigma * grad_x(2);
      }

       if (solutions.count(RADPOWER)) {
         solutions[RADPOWER][n] =
           excitonmodel->get_exciton_energy() * excitonmodel->get_real_density() /
           excitonmodel->get_radiative_recombination_rate();
       }

    }
}

// ------------- -------------- ------------
void ExcitonTransport::calc_excpol_integral(const NumericVector<Number>& x) {
  if (polaritons) {
    excpolprops.set(SimulationOptions::temperature);
    //excpolprops.calculate(x, SimulationOptions::temperature);

    do_assembly_DRSystem(1.0, NULL, NULL);
    if (excpolprops.int_g < 0) {
	    excpolprops.density_renormalization = 1;
	    excpolprops.ren_alpha = 1.0 / (1.0 - excpolprops.pol_tau * excpolprops.int_g);
	    std::cout << "G negative\n";
    } else {
    	bool success = scalarSolver.solve();
        if (!success) {
          throw ModelErrorException("Can not find solution for density renormalization."); 
	}
	excpolprops.ren_alpha = 1;
        excpolprops.density_renormalization = 1/scalarSolver.getSolution();
        do_assembly_DRSystem(scalarSolver.getSolution(), NULL, NULL); //to calculate 'int_f' last time
        std::cout << "G positive, ISOL " << scalarSolver.getSolution(); flush(std::cout);
    }
    //excpolprops.density_renormalization = 1.0 + excpolprops.int_f * excpolprops.pol_tau;
    std::cout << "excpolprops.int_f int_g= " << excpolprops.int_f << " " << excpolprops.int_g << "\n";
    std::cout << "excpolprops.density_renormalization ren_alpha= " << excpolprops.density_renormalization << " " << excpolprops.ren_alpha << "\n";

  }
}

void
ExcitonTransport::do_setup_solution_variables(void) {
  polaritons = get_option("polaritons", false);

  declare_solution_ext("generation", XGEN, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("Xpot", CHEMPOT, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("Xdens", XDENSITY, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("Xmob", XMOBILITY, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("Xsigm", XSIGMA, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("Xcurrent", J, SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "x");
  declare_solution_ext("dissociation", RDISS, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("rad_recombination", RRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("nonrad_recombination", RNONRAD, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("net_recombination", NETRECOMB, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
  declare_solution_ext("rad_power", RADPOWER, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");

  if (polaritons) {
    declare_solution_ext("n_polaritons", NPOLARITON, SolutionDescriptor::REAL, SolutionDescriptor::GLOBAL, "x");
    declare_solution_ext("exc_exc_recombination", EXCEXCPOLARITON, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
    declare_solution_ext("exc_phonon_recombination", EXCPHONONPOLARITON, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "x");
    declare_solution_ext("pol_rad_power", POLRADPOWER, SolutionDescriptor::REAL, SolutionDescriptor::GLOBAL, "W");
    declare_solution_ext("mott_density", MOTTDENSITY, SolutionDescriptor::REAL, SolutionDescriptor::GLOBAL, "x");
  }
}

void ExcitonTransport::get_solution_secure(std::map<ID, std::vector<double> >& solutions) {
  if (solutions.count(MOTTDENSITY)) {
    solutions[MOTTDENSITY][0] = excpolprops.max_density * excpolprops.lwell / (3e11); // all in cm here
  }

  if (solutions.count(NPOLARITON)) {
    std::vector<double>& solution = solutions[NPOLARITON];
    solution.resize(0);
    if (excpolprops.ren_alpha > 0.999) {
      solution.push_back(excpolprops.density_renormalization * excpolprops.int_f *excpolprops.pol_tau);
    } else {
      solution.push_back(excpolprops.ren_alpha * excpolprops.int_f *excpolprops.pol_tau);
    }
  }

  if (solutions.count(POLRADPOWER)) {
    std::vector<double>& solution = solutions[POLRADPOWER];
    solution.resize(0);

    if (excpolprops.ren_alpha > 0.999) {
      solution.push_back(excpolprops.density_renormalization * excpolprops.int_f * excpolprops.pol_energy);
    } else {
      solution.push_back(excpolprops.ren_alpha * excpolprops.int_f *excpolprops.pol_energy);
    }
  }
}

void
ExcitonTransport::do_assembly_DRSystem(const double& solution,
    double* rhs,
    double* jacobian)
{
  // We are solving 1 - I = int(g(I*exciton_density)*pol_tau)
  // 1x1 system -)
  double int_g = 0;
  double int_div_g = 0;
  double int_f = 0;

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  NumericVector<Number>& x_solution = system.get_solution_vector();

  DenseVector<Number> X;

  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  const double phi0 = get_scaling().get_potential_scaling();

  const Device& device = get_environment().get_device();
  const Options& params = get_options();

  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int ex_var = system.variable_number("fermi_x");

  FEType fe_type = system.variable_type(ex_var);

  // the finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, params.integration_order);
  fe->attach_quadrature_rule(&qrule);

  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe->get_JxW();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();
  //
  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // the DOF indices
  vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  // loop over all active elements
  for ( ; el != end_el ; ++el) {
    const Elem* elem = *el;
    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices);
    unsigned int n_dofs = dof_indices.size();

    X.resize(n_dofs);

    fe->reinit(elem);

    // extract local ddsol, accounting for constraints
    dof_map.extract_local_vector(x_solution, dof_indices, X);

    ExcitonProperties* em = dynamic_cast<ExcitonProperties*>(device.get_material(subdomain)->get_model(get_id()));

    assert(em != NULL);

    em->reinit(elem);

    //We know that Xhop is the same in element
    double Xhop = excpolprops.getXHopfield(elem, q_point[0], true);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
      // get the exciton electro-chemical potential at the quadrature point
      double ex = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++) {
        ex += phi[i][qp] * X(i);
      }

      em->set_coordinates(q_point[qp]);
      em->set_effective_potential(phi0 * ex);
      em->calculate_density();

      double x = em->get_density();

      double xx = x * solution * excpolprops.density_renormalization; // Used in f() and g(), where argument wil be divide by DR. Thats why we multiply here (!!!).

      int_g += JxW[qp] * excpolprops.g(xx, elem, q_point[qp]);
      int_div_g += JxW[qp] * excpolprops.div_g(xx, elem, q_point[qp]) * x;
      int_f += JxW[qp] * excpolprops.f(xx, elem, q_point[qp]);
    }
  }

  for (int k = 0; k < mesh.mesh_dimension(); k++) {
    int_f = int_f * get_scaling().get_length_scaling();
    int_div_g = int_div_g * get_scaling().get_length_scaling();
    int_g = int_g * get_scaling().get_length_scaling();
  }
  // We are solving 1 - I - int(g(I*exciton_density)*pol_tau) = 0;

  if (rhs != NULL) {
    *rhs = 1 - solution - excpolprops.pol_tau * int_g;
  }

  if (jacobian != NULL) {
    *jacobian = -1 - excpolprops.pol_tau * int_div_g;
  }

  if (rhs == NULL & jacobian == NULL) {
    excpolprops.int_f = int_f;
    excpolprops.int_g = int_g;
  }
}

