// $Id$

// module includes
#include "ExcitonTransport.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Scaling.h"
#include "Constants.h"
#include "ExcitonProperties.h"
#include "TiberPetscNonlinearSolver.h"

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
#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

// C++ includes

using namespace std;
using namespace DriftDiffusionDefs;


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
    mesh_units(1e-4) // default mesh units are um
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
    solver_params(rhs.solver_params),
    mesh_units(rhs.mesh_units)
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
    mesh_units = rhs.mesh_units;
  }
  return *this;
}



ExcitonTransport::SolverParameters::SolverParameters(void)
  : nonlinear_tolerance(1e-12),
    nonlinear_abs_tolerance(1e-18),
    nonlinear_step_tolerance(1e-6),
    nonlinear_max_iterations(10),
    linear_tolerance(1e-6),
    linear_abs_tolerance(1e-12),
    linear_max_iterations(500),
    ls_maxstep(0.025),
    ls_type(3),
    ksp_type(KSPBCGSL),
    pc_type(PCILU)
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

  _scaling.set_potential_scaling(phi0);
  _scaling.set_length_scaling(x0 * _options.mesh_units);
  _scaling.set_mobility_scaling(mu0);
  _scaling.set_density_scaling(C0);
}



void
ExcitonTransport::set_to_remembered_solution(void)
{
  NonlinearImplicitSystem& system =
    _eq_system->get_system<NonlinearImplicitSystem>(
        get_equation_system_name());
  *(system.solution) = system.get_vector("remembered solution");
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

  const double phi0 = _scaling.get_potential_scaling();

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
  SolverParameters& solver_params =
    get_options().solver_params;

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", 5));
  myopts.mesh_units = opts.get_option("mesh_units", 1e-4);

  myopts.mesh_refinement = opts.get_option("mesh_refinement", false);

  solver_params.nonlinear_tolerance = opts.get_option("nonlin_rel_tol", 1e-9);
  solver_params.nonlinear_abs_tolerance = opts.get_option("nonlin_abs_tol",
      1e-15);
  solver_params.nonlinear_step_tolerance = opts.get_option("nonlin_step_tol",
      1e-6);
  solver_params.nonlinear_max_iterations = opts.get_option("nonlin_max_it", 10);
  solver_params.linear_tolerance = opts.get_option("lin_rel_tol", 1e-6);
  solver_params.linear_abs_tolerance = opts.get_option("lin_abs_tol", 1e-12);
  solver_params.linear_max_iterations = opts.get_option("lin_max_it", 500);
  solver_params.ls_maxstep = opts.get_option("ls_maxstep", 0.05);
  //string pc = opts.get_option("pc", "ilu");

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

  // we can remember a solution for future use
  system.add_vector("remembered solution");
  // for adaptive mesh refinement we need the old solution
  // befor a refinement step
  system.add_vector("old solution");

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
ExcitonTransport::do_solve(void)
{

  assert(_rebuild_eq_system == false);

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;

  parse_options();

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
  NumericVector<Number>& old_solution = system.get_vector("old solution");

  // set the solver parameters (they could have change since we made
  // the first calculation)
  set_solver_params(*system.nonlinear_solver);


  try { system.solve(); }
  // if it did not converge, so try half of the step
  catch (PetscDivergedError& e)
  {
    if (e.get_solver_type() == 1) cerr << "KSP ";
    else cerr << "SNES ";
    cerr << "diverged: " << e.get_reason() <<
      " at iteration " << e.get_iteration() <<
      " (fnorm = " << e.get_fnorm() << ")\n";
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
  }

  _n_nonlinear_iterations = system.n_nonlinear_iterations();
  _final_residual = system.final_nonlinear_residual();
  old_solution = solution;



  if (params.mesh_refinement)
  {
    MeshRefinement mesh_refinement(mesh);
    ErrorVector error;
    KellyErrorEstimator error_estimator;

    bool is_refined;
    double delta_l2_norm, l2_norm;

    for (int r_step = 0; r_step < params.max_refinement_steps; r_step++)
    {
      cout << "\nRefinement step " << (r_step + 1) << "\n" << flush;
      error_estimator.estimate_error(system, error);

      mesh_refinement.flag_elements_by_error_fraction(error,
                                          params.refine_fraction,
                                          params.coarsen_fraction,
                                          params.max_refinement_level);

      is_refined = mesh_refinement.refine_and_coarsen_elements();

      equation_systems.reinit();

      if (is_refined)
        get_environment().update_element_list();


      old_solution = solution;

      try { system.solve(); }
      catch (PetscDivergedError& e)
      {
        // TODO
        if (e.get_solver_type() == 1) cerr << "KSP ";
        else cerr << "SNES ";
        cerr << "diverged: " << e.get_reason() <<
          " at iteration " << e.get_iteration() <<
          " (fnorm = " << e.get_fnorm() << ")\n";

        // we have to clear the solver context, because the exception
        // left it in a unknown state such that future solves would fail.
        solution = old_solution;
        system.nonlinear_solver->clear();
      }
      catch (PetscRuntimeError& e)
      {
        // TODO
        if (e.get_reason() == PETSC_ERR_MAT_LU_ZRPVT)
          cerr << "Zero pivot during ILU.\n";
        // we have to clear the solver context, because the exception
        // left it in a unknown state such that future solves would fail.
        solution = old_solution;
        system.nonlinear_solver->clear();
      }

/*
      solution.close();
      l2_norm = solution.l2_norm() / sqrt(mesh.n_nodes());
      old_solution.add(-1.0, solution);


      old_solution.close();
      delta_l2_norm = old_solution.l2_norm() / sqrt(mesh.n_nodes());
      cout << "  |new - old|_2 = " << delta_l2_norm
        << " |new| = " << l2_norm << "\n\n";
*/

    }

    // disable refinement
    // TODO: ???
    disable_mesh_refinement();
  }
}






double
ExcitonTransport::get_solution(const Elem* elem, const Point& p)
{
  double solution = 0.0;

  // this will contain the element in which p lies and for which
  // DriftDiffusion knows the potential
  const Elem* el = elem;

  SimulationEnvironment& env = get_environment();
  
  
  // check if elem is an active element of the simulation
  if (!env.contains_element(elem))
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();

    while ((parent != NULL) && (!env.contains_element(parent)))
      parent = parent->parent();

    el = parent; // is NULL if no parent

    // no parent, so check for children
    if (el == NULL)
    {
      vector<const Elem*> tree;
      elem->family_tree(tree, false);
      
      unsigned int len = tree.size();
      for (unsigned int i = 0; i < len; i++)
      {
        const Elem* elem_i = tree[i];
        if (env.contains_element(elem_i))
        {
          if (elem_i->contains_point(p))
          {
            // we have found it, so get out of the for loop
            el = elem_i;
            break;
          }
        }
      }
    }
  }
  // now el points to a valid element containing p or is NULL

  if (el != NULL)
    get_solution_secure(el, p, solution);

  return solution;
}





void
ExcitonTransport::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<double>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  NonlinearImplicitSystem* system;
  system = &get_equation_systems().get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

  const NumericVector<Number>& ddsol = *(system->solution);
  const NumericVector<Number>& old_ddsol = system->get_vector("old solution");
  double a = get_relaxation_factor();
  double b = 1 - a;

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("fermi_x");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

  vector<unsigned int> dof_indices_u;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_u, u_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  for (unsigned int n = 0; n < np; n++)
  {
    double u = 0;
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
      u  += phi[i][n] * (a * ddsol(dof_indices_u[i]) +
          b * old_ddsol(dof_indices_u[i]));

    // scale the potential back
    u  *= phi0;

    solution[n] = u;
  }
}






// implementation taken from libmesh equation_systems.C
void
ExcitonTransport::build_densities(vector<double>& densities,
    vector<string>& names)
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

  const unsigned int n_vars  = 4;
  names.resize(n_vars);
  names[0] = "density";
  names[1] = "net_recombination_rate";
  names[2] = "generation_rate";
  names[3] = "electro-chimical_potential";

  densities.resize(nn * n_vars);

  vector<double> local(densities.size());
  vector<unsigned short int> node_conn(nn);

  vector<double> nodal_val;

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  fill(densities.begin(), densities.end(), 0.0);
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

  MeshBase::const_element_iterator it =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_local_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

#ifdef ENABLE_INFINITE_ELEMENTS
    // infinite elements should be skipped...
    if (!elem->infinite())
#endif
    {
      dof_map.dof_indices(elem, dof_indices_u, u_var);

      ID subdomain = elem->subdomain_id();

      ExcitonProperties* excitonmodel =
        dynamic_cast<ExcitonProperties*>(
            device.get_material(subdomain)->get_model(get_id()));

      assert(excitonmodel!= NULL);

      excitonmodel->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());


      vector<Point> nodes(elem->n_nodes());
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        nodes[n] = elem->point(n);


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        double ex  = phi0 * solution(dof_indices_u[n]);
        
        excitonmodel->set_coordinates(elem->point(n));
        excitonmodel->set_effective_potential(ex);
        excitonmodel->calculate_density();
        excitonmodel->calculate_net_recombination_rate();

        assert (node_conn[elem->node(n)] != 0);

        unsigned int id = n_vars * elem->node(n);
        double nodal_val = excitonmodel->get_density();
        local[id] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = excitonmodel->get_net_recombination_rate();
        local[id + 1] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = 0;
        local[id + 2] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = ex;
        local[id + 3] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);
      }

    }
  }

#ifdef HAVE_MPI
  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  MPI_Allreduce (&local[0], &densities[0], densities.size(),
		 MPI_REAL, MPI_SUM, libMesh::COMM_WORLD);
#else
  densities = local;
#endif

}






void
ExcitonTransport::build_current_density(vector<double>& current,
    vector<string>& names)
{
/*
  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;

  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      get_equation_system_name());

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_elem();

  const unsigned int n_vars  = 6;
  names.resize(n_vars);
  names[0] = "Jn_x";
  names[1] = "Jn_y";
  names[2] = "Jn_z";
  names[3] = "Jp_x";
  names[4] = "Jp_y";
  names[5] = "Jp_z";

  current.resize(nn * n_vars);

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();
  const double x0 = get_options().mesh_units;

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, libMeshEnums::CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      device.get_element_data().get_data(elem->top_parent());
    assert(sc != NULL);

    sc->reinit(elem);

    fe->reinit(elem);

    unsigned int n_dofs = dof_indices_u.size();
    // get the solution values at the centroid
    Real en_x = 0.0, ep_x = 0.0;
    Real en_y = 0.0, ep_y = 0.0;
    Real en_z = 0.0, ep_z = 0.0;
    Real u  = 0.0;
    Real en = 0.0;
    Real ep = 0.0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      en_x  += dphi[i][0](0) * solution(dof_indices_en[i]);
      en_y  += dphi[i][0](1) * solution(dof_indices_en[i]);
      en_z  += dphi[i][0](2) * solution(dof_indices_en[i]);

      ep_x  += dphi[i][0](0) * solution(dof_indices_ep[i]);
      ep_y  += dphi[i][0](1) * solution(dof_indices_ep[i]);
      ep_z  += dphi[i][0](2) * solution(dof_indices_ep[i]);

      u  += phi[i][0] * solution(dof_indices_u[i]);
      en += phi[i][0] * solution(dof_indices_en[i]);
      ep += phi[i][0] * solution(dof_indices_ep[i]);
    }

    sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, q_point[0]);
    double sigma_e = sc->get_electron_conductivity() * Constants::e;
    double sigma_h = sc->get_hole_conductivity() * Constants::e;

    unsigned int id = n_vars * elem_number;
    //current[id] = phi0 / x0 * (sigma_e * en_x + sigma_h * ep_x);
    //current[id + 1] = phi0 / x0 * (sigma_e * en_y + sigma_h * ep_y);
    //current[id + 2] = phi0 / x0 * (sigma_e * en_z + sigma_h * ep_z);
    current[id]     = phi0 / x0 * sigma_e * en_x;
    current[id + 1] = phi0 / x0 * sigma_e * en_y;
    current[id + 2] = phi0 / x0 * sigma_e * en_z;
    current[id + 3] = phi0 / x0 * sigma_h * ep_x;
    current[id + 4] = phi0 / x0 * sigma_h * ep_y;
    current[id + 5] = phi0 / x0 * sigma_h * ep_z;

    elem_number++;
  }
  current.resize(elem_number * n_vars);
*/
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
  double x_max = 1;
  //
  // we calculate on a scaled mesh with |xmax - xmin| = 1, but we did not
  // explicitly scale the mesh, so we have to account for this in the code,
  // assuming that the mesh is drawn in units of 'mesh_units'
  const double x0_mesh = x0 / params.mesh_units;
  //
  // the scaling value for the jacobian
  double J_scale;
  switch (dim)
  {
    case 3:
      J_scale = x0_mesh * x0_mesh * x0_mesh;
      break;
    case 2:
      J_scale = x0_mesh * x0_mesh;
      break;
    default:
      J_scale = x0_mesh;
      break;
  }
  // the scaling value for the surface jacobian
  double Jface_scale = 1;
  switch (dim)
  {
    case 3:
      Jface_scale = x0_mesh * x0_mesh;
      break;
    case 2:
      Jface_scale = x0_mesh;
      break;
  }



  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int ex_var = system.variable_number("fermi_x");

  FEType fe_type = system.variable_type(ex_var);

  // the finite element
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
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

      // remember the maximum densities
      x_max = (x_max > x) ? x_max : x;

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
      double J = JxW[qp] / J_scale;

      //
      // First we will build the system matrix Ke_ij
      //
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          double laplace =
            J * (dphi[i][qp] * dphi[j][qp]) * x0_mesh * x0_mesh;

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
              double laplace = (dphi[i][qp] * dphi[k][qp]) * x0_mesh * x0_mesh;

              double elem_contrib =
                dsigma_x_phi * laplace * X(k);

              Ke(i,j) += elem_contrib;

            }

            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            Ke(i,j) += dRx * phi_i_x_phi_j ;
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

          Fe(i) += net_recomb_x;
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

