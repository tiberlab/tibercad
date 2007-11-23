
// $Id$

// module includes
#include "ExcitonTransport.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Scaling.h"
#include "Constants.h"
#include "ExcitonProperties.h"
#include "TiberPetscNonlinearSolver.h"
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
#include "nonlinear_implicit_system.h"
//#include "mesh_refinement.h"
//#include "error_vector.h"
//#include "kelly_error_estimator.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

// C++ includes

using namespace std;


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
    integration_order(rhs.integration_order),
    solver_params(rhs.solver_params)
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
    solver_params = rhs.solver_params;
  }
  return *this;
}



ExcitonTransport::SolverParameters::SolverParameters(void)
  : nonlinear_tolerance(1e-9),
    nonlinear_abs_tolerance(1e-15),
    nonlinear_step_tolerance(1e-6),
    nonlinear_max_iterations(20),
    linear_tolerance(1e-6),
    linear_abs_tolerance(1e-15),
    linear_max_iterations(500),
    ls_maxstep(0.025),
    ls_type(3),
    ksp_type(KSPBCGSL),
    pc_type(PCCOMPOSITE)
{
}

ExcitonTransport::SolverParameters::SolverParameters(const SolverParameters& rhs)
  : nonlinear_tolerance(rhs.nonlinear_tolerance),
    nonlinear_abs_tolerance(rhs.nonlinear_abs_tolerance),
    nonlinear_step_tolerance(rhs.nonlinear_step_tolerance),
    nonlinear_max_iterations(rhs.nonlinear_max_iterations),
    linear_tolerance(rhs.linear_tolerance),
    linear_abs_tolerance(rhs.linear_abs_tolerance),
    linear_max_iterations(rhs.linear_max_iterations),
    ls_maxstep(rhs.ls_maxstep),
    ls_type(rhs.ls_type),
    ksp_type(rhs.ksp_type),
    pc_type(rhs.pc_type)
{
}

ExcitonTransport::SolverParameters&
ExcitonTransport::SolverParameters::operator=(const SolverParameters& rhs)
{
  if (&rhs != this)
  {
    nonlinear_tolerance = rhs.nonlinear_tolerance;
    nonlinear_abs_tolerance = rhs.nonlinear_abs_tolerance;
    nonlinear_step_tolerance = rhs.nonlinear_step_tolerance;
    nonlinear_max_iterations = rhs.nonlinear_max_iterations;
    linear_tolerance = rhs.linear_tolerance;
    linear_abs_tolerance = rhs.linear_abs_tolerance;
    linear_max_iterations = rhs.linear_max_iterations;
    ls_maxstep = rhs.ls_maxstep;
    ls_type = rhs.ls_type;
    ksp_type = rhs.ksp_type;
    pc_type = rhs.pc_type;
  }
  return *this;
}


ExcitonTransport::ExcitonTransport(void)
  : _rebuild_eq_system(true)
{
}


ExcitonTransport::~ExcitonTransport(void)
{
  cleanup_solver();
}




PhysicalModel*
ExcitonTransport::create_physical_model(const ModelOptions& options) const
throw (ModelErrorException)
{
  const string& modelname = options.get_option("model", "simple");

  ExcitonProperties* model =
    ExcitonProperties::create(modelname, options);

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

  const Mesh& mesh = get_mesh();
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
    get_equation_systems().delete_system(get_equation_system_name());
    _rebuild_eq_system = true;
  }
}


void
ExcitonTransport::cleanup_solver(void)
{
  reset_solver();
}


void
ExcitonTransport::set_solver_params(NonlinearSolver<Number>& solver)
{
  SolverClass& solver_class =
    static_cast<SolverClass&>(solver);

  SolverParameters& solver_params = _options.solver_params;

  const double phi0 = get_scaling().get_potential_scaling();

  unsigned int nonlin_max_its = solver_params.nonlinear_max_iterations;

  double sqrt_n = sqrt((double) get_mesh().n_nodes());

  solver_class.set_snes_options(solver_params.nonlinear_tolerance,
      solver_params.nonlinear_abs_tolerance * sqrt_n,
      solver_params.nonlinear_step_tolerance, nonlin_max_its);
  
  solver_class.set_snes_ls_options(solver_params.ls_type,
      solver_params.ls_maxstep * sqrt_n / phi0);

  solver_class.set_ksp_options(solver_params.linear_tolerance,
      solver_params.linear_abs_tolerance * sqrt_n,
      solver_params.linear_max_iterations);

  solver_class.set_ksp_type(solver_params.ksp_type);
  solver_class.set_pc_type(solver_params.pc_type);
}




void
ExcitonTransport::parse_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();
  SolverParameters& solver_params = myopts.solver_params;

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", (int) myopts.integration_order));

  myopts.mesh_refinement = opts.get_option("mesh_refinement", false);

  solver_params.nonlinear_tolerance = opts.get_option("nonlin_rel_tol",
      solver_params.nonlinear_tolerance);
  solver_params.nonlinear_abs_tolerance = opts.get_option("nonlin_abs_tol",
      solver_params.nonlinear_abs_tolerance);
  solver_params.nonlinear_step_tolerance = opts.get_option("nonlin_step_tol",
      solver_params.nonlinear_step_tolerance);
  solver_params.nonlinear_max_iterations = opts.get_option("nonlin_max_it",
      solver_params.nonlinear_max_iterations);
  solver_params.linear_tolerance = opts.get_option("lin_rel_tol",
      solver_params.linear_tolerance);
  solver_params.linear_abs_tolerance = opts.get_option("lin_abs_tol",
      solver_params.linear_abs_tolerance);
  solver_params.linear_max_iterations = opts.get_option("lin_max_it",
      solver_params.linear_max_iterations);
  solver_params.ls_maxstep = opts.get_option("ls_maxstep",
      solver_params.ls_maxstep);

  string pc = opts.get_option("pc_type", "");
  if (pc == "") {}
  else if (pc == "ilu")
    solver_params.pc_type = PCILU;
  else if (pc == "composite")
    solver_params.pc_type = PCCOMPOSITE;
  else if (pc == "jacobi")
    solver_params.pc_type = PCJACOBI;
  else if (pc == "lu")
    solver_params.pc_type = PCLU;
  else if (pc == "cholesky")
    solver_params.pc_type = PCCHOLESKY;


}

void
ExcitonTransport::do_init(void)
{
  if (!_rebuild_eq_system) return;

  _device = &get_environment().get_device();

  EquationSystems& equation_systems = get_equation_systems();

  // create the exciton continuity equation
  NonlinearImplicitSystem& system =
    equation_systems.add_system<NonlinearImplicitSystem>(
        get_equation_system_name());


  // we use PETSc
  system.nonlinear_solver =
    AutoPtr<NonlinearSolver<Number> >(new SolverClass);

  // set the options for the PETSc nonlinear solver
  set_solver_params(*system.nonlinear_solver);

  system.add_variable("fermi_x", libMeshEnums::FIRST);

  // a local scaling
  system.add_vector("local_scaling");

  system.nonlinear_solver->matvec = assemble;

  // set some parameters (but we don't use them in this way)
  // they have to exist for libmesh 
  SolverParameters& solver_params = get_options().solver_params;

  equation_systems.parameters.set<unsigned int>(
    "nonlinear solver maximum iterations") = 100;
      
  equation_systems.parameters.set<Real>("nonlinear solver tolerance") = 1e-6;


  // initialize the newly created system
  system.init();

  // compute scaling factors
  compute_scaling();


  _rebuild_eq_system = false;

}





void
ExcitonTransport::set_initial_guess(double guess)
{
  assert(_rebuild_eq_system == false);

  NonlinearImplicitSystem& system =
    get_equation_systems().get_system<NonlinearImplicitSystem>(
        get_equation_system_name());

  NumericVector<Number>& solution = *(system.solution);

  solution.zero();
  solution.add(guess);
}



void
ExcitonTransport::build_local_scaling(void)
{
  NonlinearImplicitSystem* system = 
    &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());
  //TiberNonlinearSystem* system =
  //  &get_equation_systems().get_system<TiberNonlinearSystem>(
  //      get_equation_system_name());
  
  
  const NumericVector<Number>& solution = *(system->solution);
  //const NumericVector<Number>& solution = system->get_solution_vector();
  NumericVector<Number>& locscal = system->get_vector("local_scaling");

  // aliases for nicer code
  const Device& device = *(_device);
  const Mesh& mesh = get_mesh();

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

    //assert(elem->n_nodes() == dof_indices_u.size());


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


      double sigma_x = JxW[qp] * sc->get_mobility() * sc->get_density();

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //if (sigma_x > 1)
        //  locscal.add(dof_indices[i], sigma_x * (dphi[i][qp] * dphi[i][qp]));
        //else
          locscal.set(dof_indices[i], 1.0);
      }


    } // end loop over quadrature points
  } // end loop over elements

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

  // aliases for nicer code
  Options& params = get_options();
  SolverParameters& solver_params = params.solver_params;

  Mesh& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  EquationSystems& equation_systems = get_equation_systems();

  if (dim == 1)
    if (solver_params.ksp_type == KSPBCGSL)
      solver_params.ksp_type = KSPBCGS;


  NonlinearImplicitSystem& system =
    equation_systems.get_system<NonlinearImplicitSystem>(
        get_equation_system_name());

  NumericVector<Number>& solution = *(system.solution);

  // set the solver parameters (they could have change since we made
  // the first calculation)
  set_solver_params(*system.nonlinear_solver);


  bool failure = true;
  string msg("ExcitonTransport: solve failed (");

  try
  {
    system.solve();

    failure = false;
  }
  catch (PetscDivergedError& e)
  {
    if (e.get_solver_type() == 1) cerr << "KSP ";
    else cerr << "SNES ";
    cerr << "diverged: " << e.get_reason() <<
      " at iteration " << e.get_iteration() <<
      " (fnorm = " << e.get_fnorm() << ")\n";

    msg += e.what();
    msg += ")\n";
  }
  catch (PetscRuntimeError& e)
  {
    cerr << "Petsc runtime error: " << e.get_reason();
    if (e.get_reason() == PETSC_ERR_MAT_LU_ZRPVT)
    {
      // in the case of a zero pivot in (I)LU factorization
      // we try another preconditioner
      cerr << " (Zero pivot during ILU.)";
      PCType old_pc = solver_params.pc_type;
      solver_params.pc_type = PCJACOBI;
      set_solver_params(*system.nonlinear_solver);
      solver_params.pc_type = old_pc;
    }
    cerr << "\n";
    msg += e.what();
    msg += ")\n";
  }

  _n_nonlinear_iterations = system.n_nonlinear_iterations();
  _final_residual = system.final_nonlinear_residual();

  if (failure)
  {
    system.nonlinear_solver->clear();
    do_init();
    throw SolveFailedException(msg);
  }

}





// implementation taken from libmesh equation_systems.C
void
ExcitonTransport::build_nodal_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{
  
  NonlinearImplicitSystem* system;

  system = &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

  // aliases for nicer code
  const Device& device = *(_device);
  const Mesh& mesh = get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();
  
  legend.resize(variables.size());

  // for each possible variable we set the vector index
  // -1 means, the variable should not be plotted
  unsigned int n_vars = 0;
  const set<string>::const_iterator varend(variables.end());


  int Ef = -1;
  if (variables.find("xEffPot") != varend)
  {
    Ef = n_vars;
    legend[n_vars] = "exciton_effective_potential";
    n_vars++;
  }


  int xdens = -1;
  if (variables.find("xDensity") != varend)
  {
    xdens = n_vars;
    legend[n_vars] = "exciton_density";
    n_vars++;
  }


  int rec = -1;
  int num_rec = 0;
  map<ID, string> rec_model_ids;
  if (variables.find("ExcitonRecombination") != varend)
  {
    num_rec = 5;
    rec = n_vars;
    
    legend.resize(variables.size() + num_rec);

    legend[n_vars] = "exciton_generation";
    n_vars++;

    legend[n_vars] = "exciton_radiative_recombination";
    n_vars++;

    legend[n_vars] = "exciton_nonradiative_recombination";
    n_vars++;

    legend[n_vars] = "exciton_dissociation";
    n_vars++;

    legend[n_vars] = "net_exciton_recombination";
    n_vars++;
  }


  int mu = -1;
  if (variables.find("xMob") != varend)
  {
    mu = n_vars;
    legend[n_vars] = "exciton_mobility";
    n_vars++;
  }

  legend.resize(n_vars);
    
  results.resize(nn * n_vars);

  vector<double> local(results.size());
  vector<unsigned short int> node_conn(nn);

  vector<double> nodal_val;

  // the scaling parameters to scale back the result
  const double phi0 = get_scaling().get_potential_scaling();


  fill(results.begin(), results.end(), 0.0);
  fill(local.begin(), local.end(), 0.0);

  // Get the number of elements that share each node.  We will
  // compute the average value at each node.
  {
    vector<unsigned short int> node_conn_local(node_conn.size());
    
    MeshBase::const_element_iterator it =
      mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end =
      mesh.active_local_elements_end(); 

    for ( ; it != end; ++it)
      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	node_conn_local[(*it)->node(n)]++;

#ifdef HAVE_MPI
    // Gather the distributed node_conn arrays in the case of
    // multiple processors
    //
    // (Note that we use an unsigned short int here even though an
    // unsigned char would be more that sufficient.  The MPI 1.1
    // standard does not require that MPI_SUM, MPI_PROD etc... be
    // implemented for char data types. 12/23/2003 - BSK)  
    MPI_Allreduce (&node_conn_local[0], &node_conn[0], node_conn.size(),
		   MPI_UNSIGNED_SHORT, MPI_SUM, libMesh::COMM_WORLD);
    
#else
    // Without MPI the node_conn_local and the node_conn arrays
    // are necessarily identical
    node_conn = node_conn_local;
    
#endif
  }

  const unsigned int u_var = system->variable_number("fermi_x");
  
  vector<unsigned int> dof_indices_u;

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  QGauss qrule(dim, libMeshEnums::CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  MeshBase::const_element_iterator it =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_local_elements_end(); 

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

#ifdef ENABLE_INFINITE_ELEMENTS
    // infinite elements should be skipped...
    if (!elem->infinite())
#endif
    {
      dof_map.dof_indices(elem, dof_indices_u, u_var);

      ExcitonProperties* excitonmodel =
        dynamic_cast<ExcitonProperties*>(
            device.get_material(subdomain)->get_model(get_id()));

      assert(excitonmodel!= NULL);

      excitonmodel->reinit(elem);

      fe->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        double ex  = phi0 * solution(dof_indices_u[n]);

        double kT = excitonmodel->get_lattice_temperature();
        
        excitonmodel->set_carrier_temperature(kT);
        excitonmodel->set_coordinates(elem->point(n));
        excitonmodel->set_effective_potential(ex);
        excitonmodel->calculate_density();
        excitonmodel->calculate_net_recombination_rate();


        assert (node_conn[elem->node(n)] != 0);
        double conn = static_cast<double>(node_conn[elem->node(n)]);

        unsigned int id = n_vars * elem->node(n);

        if (xdens != -1)
        {
          double nodal_val = excitonmodel->get_density();
          local[id + xdens] += nodal_val / conn;
        }

        if (rec != -1)
        {
          // recombination models
          int first = id + rec;
          double tot = excitonmodel->get_net_recombination_rate();
          local[first + 4] += tot / conn;
          
          double nodal_val = excitonmodel->get_radiative_recombination_rate();
          local[first + 1] += nodal_val / conn;
          double totrec = nodal_val;
          
          nodal_val = excitonmodel->get_nonradiative_recombination_rate();
          local[first + 2] += nodal_val / conn;
          totrec += nodal_val;
          
          nodal_val = excitonmodel->get_dissociation_rate();
          local[first + 3] += nodal_val / conn;
          totrec += nodal_val;

          // note the sign: we want to plot it as positive quantity
          nodal_val = (totrec - tot) / conn;
          local[first] += nodal_val; // generation
        }


        if (mu != -1)
        {
          double nodal_val = excitonmodel->get_mobility();
          local[id + mu] += nodal_val / conn;
        }

        if (Ef != -1)
          local[id + Ef] += ex / conn;


      } // end loop over nodes
    } // end if not infinite element
  } // end loop over elements

#ifdef HAVE_MPI
  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  MPI_Allreduce (&local[0], &results[0], results.size(),
		 MPI_REAL, MPI_SUM, libMesh::COMM_WORLD);
#else
  results = local;
#endif

}






void
ExcitonTransport::build_elemental_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{

  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;
  
  NonlinearImplicitSystem* system;

  system = &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

  // aliases for nicer code
  const Device& device = *(_device);
  const Mesh& mesh = get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_elem();
  
  legend.resize(variables.size());

  // for each possible variable we set the vector index
  // -1 means, the variable should not be plotted
  unsigned int n_vars = 0;
  const set<string>::const_iterator varend(variables.end());


  int J = -1;
  if (variables.find("xCurrent") != varend)
  {
    J = n_vars;
    legend.resize(variables.size() + dim);
    switch (dim)
    {
      case 3:
        legend[J + 2] = "Jx_z[cm^2*s^-1]";
        n_vars++;
      case 2:
        legend[J + 1] = "Jx_y[cm^2*s^-1]";
        n_vars++;
        legend[J + dim] = "|J|[cm^2*s^-1]";
        n_vars++;
      default:
        legend[J] = "Jx_x[cm^-2*s^-1)]";
        n_vars++;
    }
  }


  legend.resize(n_vars);

  results.resize(nn * n_vars);

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  const unsigned int u_var = system->variable_number("fermi_x");
  
  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  QGauss qrule(dim, libMeshEnums::CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  vector<unsigned int> dof_indices_u;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    ExcitonProperties* excitonmodel =
      dynamic_cast<ExcitonProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(excitonmodel!= NULL);

    excitonmodel->reinit(elem);

    fe->reinit(elem);
    
    unsigned int n_dofs = dof_indices_u.size();
    // get the solution values at the centroid
    Real ex_x = 0.0;
    Real ex_y = 0.0;
    Real ex_z = 0.0;
    Real ex = 0.0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      ex_x  += dphi[i][0](0) * solution(dof_indices_u[i]);
      ex_y  += dphi[i][0](1) * solution(dof_indices_u[i]);
      ex_z  += dphi[i][0](2) * solution(dof_indices_u[i]);
      
      ex += phi[i][0] * solution(dof_indices_u[i]);
    }
    ex_x *= phi0;
    ex_y *= phi0;
    ex_z *= phi0;

    // prepare for calculating local properties
    excitonmodel->set_coordinates(elem->centroid());

    double kT = excitonmodel->get_lattice_temperature();

    excitonmodel->set_carrier_temperature(kT);
    excitonmodel->set_coordinates(elem->centroid());
    excitonmodel->set_effective_potential(phi0 * ex);
    excitonmodel->calculate_density();


    double sigma = excitonmodel->get_density() * excitonmodel->get_mobility();

    unsigned int id = n_vars * elem_number;


    if (J != -1)
    {
      double jx = sigma * ex_x;
      double jy = sigma * ex_y;
      double jz = sigma * ex_z;
      switch (dim)
      {
        case 3:
          results[id + J + 2] = jz;
        case 2:
          results[id + J + 1] = jy;
          results[id + J + dim] = sqrt(jx * jx + jy * jy + jz * jz);
        default:
          results[id + J] = jx;
      }
    }


    elem_number++;
  }

  results.resize(elem_number * n_vars);
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
  const Mesh& mesh = get_mesh();

  EquationSystems& eq_sys = get_equation_systems();
  NonlinearImplicitSystem& system =
    eq_sys.get_system<NonlinearImplicitSystem>(get_equation_system_name());

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
      double sigma_x = (mux * x) / (mu0 * C0);


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

            Ke(i,j) += sigma_x * laplace ;

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
}





ID
ExcitonTransport::convert_variable_name_to_id(const string& variable_name) const
{
  ID id = INVALID_ID;


  if (variable_name == "chemPot")
    id = CHEMPOT;
  else if (variable_name == "xCond")
    id = XSIGMA;
  else if (variable_name == "xMob")
    id = XMOBILITY;
  else if (variable_name == "xDensity")
    id = XDENSITY;
  else if (variable_name == "J")
    id = J;
  else if (variable_name == "J_x")
    id = JX;
  else if (variable_name == "J_y")
    id = JY;
  else if (variable_name == "J_z")
    id = JZ;
  else if (variable_name == "Rad_power")
    id = RADPOWER;

  else if (variable_name == "dissociation")
    id = RDISS;



  return id;
}

      


void
ExcitonTransport::get_solution_secure(const Elem* elem,
    const set<ID>& ids, vector<map<ID, double> >& values)
{
  unsigned int np = elem->n_nodes();

  vector<Point> points(np);
  for (int i = 0; i < np; i++)
    points[i] = elem->local_node(elem->type(), i);


  NonlinearImplicitSystem* system;
  system = &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

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

  for (unsigned int n = 0; n < np; n++)
  {
    // the chemical potential
    double ex = phi0 * ddsol(dof_indices_u[n]);

    // do interpolation
    double grad_x = 0.0;
    double grad_y = 0.0;
    double grad_z = 0.0;
    double dT_x = 0.0;
    double dT_y = 0.0;
    double dT_z = 0.0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      grad_x  += dphi[i][n](0) * ddsol(dof_indices_u[i]);
      grad_y  += dphi[i][n](1) * ddsol(dof_indices_u[i]);
      grad_z  += dphi[i][n](2) * ddsol(dof_indices_u[i]);
      
      //dT_x  += dphi[i][n](0) * T_nodes[i];
      //dT_y  += dphi[i][n](1) * T_nodes[i];
      //dT_z  += dphi[i][n](2) * T_nodes[i];
    }

    // scale the potential back
    grad_x *= phi0;
    grad_y *= phi0;
    grad_z *= phi0;


    // prepare for calculating local properties
    excitonmodel->set_coordinates(points[n]);

    double T_lat = excitonmodel->get_lattice_temperature();
    // all are at lattice temperature
    excitonmodel->set_carrier_temperature(T_lat);

    excitonmodel->set_effective_potential(ex);

    excitonmodel->calculate_density();

    double sigma = excitonmodel->get_density() * excitonmodel->get_mobility();


    if (ids.count(CHEMPOT))
      values[n][CHEMPOT] = ex;

    if (ids.count(XSIGMA))
      values[n][XSIGMA] = sigma;

    if (ids.count(XMOBILITY))
      values[n][XMOBILITY] = excitonmodel->get_mobility();

    if (ids.count(XDENSITY))
      values[n][XDENSITY] = excitonmodel->get_density();

    if (ids.count(RDISS))
      values[n][RDISS] = excitonmodel->get_dissociation_rate();

    if (ids.count(J))
    {
      double tmp = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
      values[n][J] = sigma * sqrt(tmp);
    }

    if (ids.count(JX))
      values[n][JX] = sigma * grad_x;

    if (ids.count(JY))
      values[n][JY] = sigma * grad_y;

    if (ids.count(JZ))
      values[n][JZ] = sigma * grad_z;

     if (ids.count(RADPOWER))
       values[n][RADPOWER] =  excitonmodel->get_density()/excitonmodel->get_radiative_recombination_rate() *
	 (excitonmodel->get_exciton_energy());

  }

}

      


void
ExcitonTransport::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{
  unsigned int np = p.size();

  NonlinearImplicitSystem* system;
  system = &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

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

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);

  fe->reinit(elem, &points);

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

  for (unsigned int n = 0; n < np; n++)
  {
    
    // do interpolation
    double ex = 0;
    double grad_x = 0.0;
    double grad_y = 0.0;
    double grad_z = 0.0;
    double dT_x = 0.0;
    double dT_y = 0.0;
    double dT_z = 0.0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      ex  += phi[i][n] * ddsol(dof_indices_u[i]);

      grad_x  += dphi[i][n](0) * ddsol(dof_indices_u[i]);
      grad_y  += dphi[i][n](1) * ddsol(dof_indices_u[i]);
      grad_z  += dphi[i][n](2) * ddsol(dof_indices_u[i]);
      
      //dT_x  += dphi[i][n](0) * T_nodes[i];
      //dT_y  += dphi[i][n](1) * T_nodes[i];
      //dT_z  += dphi[i][n](2) * T_nodes[i];
    }

    // scale the potential back
    ex  *= phi0;
    grad_x *= phi0;
    grad_y *= phi0;
    grad_z *= phi0;


    // prepare for calculating local properties
    excitonmodel->set_coordinates(points[n]);

    double T_lat = excitonmodel->get_lattice_temperature();
    // all are at lattice temperature
    excitonmodel->set_carrier_temperature(T_lat);

    excitonmodel->set_effective_potential(ex);

    excitonmodel->calculate_density();

    double sigma = excitonmodel->get_density() * excitonmodel->get_mobility();


    if (ids.count(CHEMPOT))
      values[n][CHEMPOT] = ex;

    if (ids.count(XSIGMA))
      values[n][XSIGMA] = sigma;

    if (ids.count(XMOBILITY))
      values[n][XMOBILITY] = excitonmodel->get_mobility();

    if (ids.count(XDENSITY))
      values[n][XDENSITY] = excitonmodel->get_density();

    if (ids.count(RDISS))
      values[n][RDISS] = excitonmodel->get_dissociation_rate();

    if (ids.count(J))
    {
      double tmp = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
      values[n][J] = sigma * sqrt(tmp);
    }

    if (ids.count(JX))
      values[n][JX] = sigma * grad_x;

    if (ids.count(JY))
      values[n][JY] = sigma * grad_y;

    if (ids.count(JZ))
      values[n][JZ] = sigma * grad_z;
  }

}



