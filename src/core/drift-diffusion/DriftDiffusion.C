// $Id$

// module includes
#include "DriftDiffusion.h"
#include "DDevice.h"
#include "Scaling.h"
#include "ElementData.h"
#include "BoundaryData.h"
#include "BoundaryDescriptor.h"
#include "Constants.h"
#include "SemiconductorModel.h"
#include "TiberPetscNonlinearSolver.h"

// libmesh includes
#include "node.h"
#include "mesh.h"
#include "mesh_modification.h"
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
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "gmv_io.h"

// C++ includes

using namespace std;


DriftDiffusion*
DriftDiffusion::_this;


DriftDiffusion::Options::Options(void)
  : mesh_refinement(false),
    max_refinement_steps(5),
    max_refinement_level(8),
    refine_fraction(0.7),
    coarsen_fraction(0.3),
    refinement_tolerance(1e-6),
    min_voltage_step(2e-3),
    integration_order(libMeshEnums::FIFTH),
    approximation_order(libMeshEnums::FIRST),
    solver_method(NEWTON),
    max_gummel_iterations(5),
    mesh_units(1e-4), // default mesh units are um
    scaling_type(Scaling::UNITS),
    n_max(0), p_max(0),
    C0_e(1), C0_h(1)
{
}

DriftDiffusion::Options::Options(const Options& rhs)
  : mesh_refinement(rhs.mesh_refinement),
    max_refinement_steps(rhs.max_refinement_steps),
    max_refinement_level(rhs.max_refinement_level),
    refine_fraction(rhs.refine_fraction),
    coarsen_fraction(rhs.coarsen_fraction),
    refinement_tolerance(rhs.refinement_tolerance),
    min_voltage_step(rhs.min_voltage_step),
    integration_order(rhs.integration_order),
    approximation_order(rhs.approximation_order),
    solver_method(rhs.solver_method),
    max_gummel_iterations(rhs.max_gummel_iterations),
    solver_params(rhs.solver_params),
    mesh_units(rhs.mesh_units),
    scaling_type(rhs.scaling_type),
    n_max(rhs.n_max),
    p_max(rhs.p_max),
    C0_e(rhs.C0_e),
    C0_h(rhs.C0_h)
{
}

DriftDiffusion::Options&
DriftDiffusion::Options::operator=(const Options& rhs)
{
  if (&rhs != this)
  {
    mesh_refinement = rhs.mesh_refinement;
    max_refinement_steps = rhs.max_refinement_steps;
    max_refinement_level = rhs.max_refinement_level;
    refine_fraction = rhs.refine_fraction;
    coarsen_fraction = rhs.coarsen_fraction;
    refinement_tolerance = rhs.refinement_tolerance;
    min_voltage_step = rhs.min_voltage_step;
    integration_order = rhs.integration_order;
    approximation_order = rhs.approximation_order;
    solver_method = rhs.solver_method;
    max_gummel_iterations = rhs.max_gummel_iterations;
    solver_params = rhs.solver_params;
    mesh_units = rhs.mesh_units;
    scaling_type = rhs.scaling_type;
    n_max = rhs.n_max;
    p_max = rhs.p_max;
    C0_e = rhs.C0_e;
    C0_h = rhs.C0_h;
  }
  return *this;
}



DriftDiffusion::SolverParameters::SolverParameters(void)
  : nonlinear_tolerance(1e-12), 
    nonlinear_abs_tolerance(1e-18),
    nonlinear_max_iterations(10),
    linear_tolerance(1e-6),
    linear_abs_tolerance(1e-12),
    linear_max_iterations(500),
    ls_maxstep(0.025),
    ksp_type(KSPBCGSL),
    pc_type(PCJACOBI)
{
}

DriftDiffusion::SolverParameters::SolverParameters(const SolverParameters& rhs)
  : nonlinear_tolerance(rhs.nonlinear_tolerance), 
    nonlinear_abs_tolerance(rhs.nonlinear_abs_tolerance),
    nonlinear_max_iterations(rhs.nonlinear_max_iterations),
    linear_tolerance(rhs.linear_tolerance),
    linear_abs_tolerance(rhs.linear_abs_tolerance),
    linear_max_iterations(rhs.linear_max_iterations),
    ls_maxstep(rhs.ls_maxstep),
    ksp_type(rhs.ksp_type),
    pc_type(rhs.pc_type)
{
}

DriftDiffusion::SolverParameters&
DriftDiffusion::SolverParameters::operator=(const SolverParameters& rhs)
{
  if (&rhs != this)
  {
    nonlinear_tolerance = rhs.nonlinear_tolerance; 
    nonlinear_abs_tolerance = rhs.nonlinear_abs_tolerance;
    nonlinear_max_iterations = rhs.nonlinear_max_iterations;
    linear_tolerance = rhs.linear_tolerance;
    linear_abs_tolerance = rhs.linear_abs_tolerance;
    linear_max_iterations = rhs.linear_max_iterations;
    ls_maxstep = rhs.ls_maxstep;
    ksp_type = rhs.ksp_type;
    pc_type = rhs.pc_type;
  }
  return *this;
}



// TODO add integrity test
DriftDiffusion::DriftDiffusion(DD::Device* device)
  : _eq_system(NULL),
    _rebuild_eq_system(true)
{
  // should throw exception
  assert(device->check_integrity());
  _device = device;

  find_dirichlet_nodes();
}

// TODO add integrity test
DriftDiffusion::DriftDiffusion(DD::Device* device,
    DriftDiffusion::Options& params)
  : _eq_system(NULL),
    _rebuild_eq_system(true)
{
  // should throw exception
  assert(device->check_integrity());
  _device = device;
  _options = params;

  find_dirichlet_nodes();
}

DriftDiffusion::~DriftDiffusion(void)
{
  delete _eq_system;
}
 
// TODO add integrity test
void
DriftDiffusion::set_device(DD::Device* device)
{
  // should throw exception
  assert(device->check_integrity());

  cleanup_solver();
    
  _device = device;

  find_dirichlet_nodes();

}


void
DriftDiffusion::compute_scaling(Scaling::ScalingType type)
{
  if (type == Scaling::NONE)
  {
    _scaling.set_scaling_type(type);
    _scaling.set_potential_scaling(1);
    _scaling.set_length_scaling(1);
    _scaling.set_mobility_scaling(1);
    _scaling.set_density_scaling(1);
    return;
  }
  
  // the scaling parameters should never be zero
  // they are in any case positive, so it will
  // always find the maximum
  double x0 = -1;
  double phi0 = -1;
  double mu0 = -1;
  double C0 = -1;
  double ni0 = 1; // let 1 be the minimum for ni0
  double eps0 = -1;
  
  // find minimum or maximum by looping over all elements
  set<SemiconductorModel*>& materials = _device->get_materials();
  DD::Device::const_material_iterator it = materials.begin();
  const DD::Device::const_material_iterator end = materials.end();
  while (it != end)
  {
    SemiconductorModel* sc = *it;
    const SemiconductorModel::MaterialDescriptor& material =
      sc->get_material_descriptor();
    
    double vt = sc->get_thermal_voltage();
    phi0 = (phi0 > vt) ? phi0 : vt;
    
    double mu = material.hole_mobility;
    mu0 = (mu0 > mu) ? mu0 : mu;
    mu = material.electron_mobility;
    mu0 = (mu0 > mu) ? mu0 : mu;

    double C = fabs(material.n_dopant.get_doping_density() -
        material.p_dopant.get_doping_density());
    C0 = (C0 > C) ? C0 : C;

    double ni = sc->get_intrinsic_density();
    ni0 = (ni0 > ni) ? ni0 : ni;

    double eps = material.relative_permittivity;
    eps0 = (eps0 > eps) ? eps0 : eps;

    ++it;
  }

  switch (type)
  {
    case Scaling::DEMARI:
      C0 = ni0;
      mu0 = 1.0 / phi0;
      // the factor 1e-2 is here because we calculate in cm
      x0 = sqrt(eps0 * 1e-2 * phi0 / (Constants::e * ni0));
      break;
      
    default: // UNITS
      C0 = (C0 > ni0) ? C0 : ni0;

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

      break;
  }

  _scaling.set_scaling_type(type);
  _scaling.set_potential_scaling(phi0);
  _scaling.set_length_scaling(x0 * _options.mesh_units);
  _scaling.set_mobility_scaling(mu0);
  _scaling.set_density_scaling(C0);
}


void
DriftDiffusion::set_simulation_voltage(const string& boundary,
    double voltage)
{
  BoundaryDescriptor* desc = _device->get_boundary(boundary);

  if (desc != NULL)
  {
    double phi0 = get_scaling().get_potential_scaling();

    _simulation_voltages[desc] = voltage;
  }
}

void
DriftDiffusion::remember_current_solution(void)
{
  _remembered_voltages = _simulation_voltages;

  NonlinearImplicitSystem& system =
    _eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
  system.get_vector("remembered solution") = *(system.solution);
  system.get_vector("remembered solution").close();
}

void
DriftDiffusion::set_to_remembered_solution(void)
{
  assert(_old_sim_voltages.size() == _remembered_voltages.size());
  _old_sim_voltages = _remembered_voltages;

  NonlinearImplicitSystem& system =
    _eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
  *(system.solution) = system.get_vector("remembered solution");

}

int
DriftDiffusion::find_boundary(const Elem* elem, int side,
    const Elem* top_parent)
{
  int top_side = -1;

  const DD::Device& device = _this->get_device();

  vector<int> top_sides =
    device.get_boundary_data().find_element(top_parent);

  if (top_sides.size() != 0)
  {
    // side numbers are inherited from the parent element, so we
    // can use them for the check
    vector<int>::iterator end = top_sides.end();
    if (find(top_sides.begin(), end, side) != end)
      top_side = side;
  }

  return top_side;
}


void
DriftDiffusion::find_dirichlet_nodes(void)
{

  const BoundaryData& bound_data = _device->get_boundary_data();
  int dim = _device->get_mesh().mesh_dimension();

  // we make a list which contains only dirichlet type boundaries
  // this should make the rest somewhat faster
  set<BoundaryDescriptor*> dirichlet_boundaries;

  DD::Device::BoundaryList& boundaries = _device->get_boundaries();
  DD::Device::const_boundary_iterator it = boundaries.begin();
  const DD::Device::const_boundary_iterator end = boundaries.end();
  for ( ; it != end; ++it)
  {
    if (((*it)->get_type("potential") == BoundaryDescriptor::DIRICHLET)
        || ((*it)->get_type("fermi_e") == BoundaryDescriptor::DIRICHLET)
        || ((*it)->get_type("fermi_h") == BoundaryDescriptor::DIRICHLET))
      dirichlet_boundaries.insert(*it);
  }

  set<BoundaryDescriptor*>::iterator not_dirichlet =
    dirichlet_boundaries.end();

  // loop over all level 0 boundary elements
  BoundaryData::const_iterator bd_it = bound_data.sides_begin();
  const BoundaryData::const_iterator bd_end = bound_data.sides_end();
  for ( ; bd_it != bd_end; ++bd_it)
  {
    if (dirichlet_boundaries.find(bd_it->second) != not_dirichlet)
    {
      const Elem* elem = (bd_it->first).first;
      int side_num = (bd_it->first).second;

      // for 1D we are done as no new boundary nodes can be added
      // side_num in this case is the node number
      if (dim == 1)
      {
        (_dirichlet_nodes)[elem->get_node(side_num)] = bd_it->second;

        continue;
      }
      
      // get the active family tree of this element
      vector<const Elem*> fam_tree;
      elem->active_family_tree(fam_tree);

      // loop over all active children and find boundary sides that
      // correspond to side_num
      vector<const Elem*>::const_iterator elem_it;
      for (elem_it = fam_tree.begin();
           elem_it != fam_tree.end(); ++elem_it)
      {
        const Elem* child = *elem_it;

        if (child->neighbor(side_num) == NULL)
        {
            AutoPtr<Elem> side = child->build_side(side_num);
            for (int i = 0; i < side->n_nodes(); i++)
              (_dirichlet_nodes)[side->get_node(i)] = bd_it->second;
        }
      }
    }
  }
}

void
DriftDiffusion::prepare_solver(void)
{

  if (_rebuild_eq_system)
  {
    // we assume that _eq_system was deleted before
    assert(_eq_system == NULL);
    _eq_system = new EquationSystems(_device->get_mesh());

    DD::Device::BoundaryList& boundaries = _device->get_boundaries();
    DD::Device::const_boundary_iterator it = boundaries.begin();
    const DD::Device::const_boundary_iterator end = boundaries.end();
    const map<const BoundaryDescriptor*,
          double>::const_iterator bc_end = _simulation_voltages.end();
    for ( ; it != end; ++it)
    {
      if (_simulation_voltages.find(*it) == bc_end)
        _simulation_voltages[*it] = 0.0;

      _old_sim_voltages[*it] = 0.0;
    }
  }

  _this = this;

  compute_scaling(_options.scaling_type);
}

void
DriftDiffusion::initialize_eq_system(EquationSystems& system)
{

  SolverParameters& solver_params =
    (_this->get_options()).solver_params;

  system.parameters.set<unsigned int>(
    "nonlinear solver maximum iterations") = 
      solver_params.nonlinear_max_iterations;
  
  system.parameters.set<Real>("nonlinear solver tolerance") =
    solver_params.nonlinear_tolerance;

  system.parameters.set<Real>("linear solver tolerance") =
    solver_params.linear_tolerance;

  system.init();

}

void
DriftDiffusion::reset_solver(void)
{
  delete _eq_system;
  _eq_system = NULL;

  _rebuild_eq_system = true;
}

void
DriftDiffusion::cleanup_solver(void)
{
  delete _eq_system;
  _eq_system = NULL;

  // erase old simulation voltage data structure
  _old_sim_voltages.erase(_old_sim_voltages.begin(),
      _old_sim_voltages.end());

  // erase simulation voltage data structure
  _simulation_voltages.erase(_simulation_voltages.begin(),
      _simulation_voltages.end());

  // erase remembered voltage data structure
  _remembered_voltages.erase(_remembered_voltages.begin(),
      _remembered_voltages.end());

  // erase boundary current data structure
  _boundary_currents.erase(_boundary_currents.begin(),
      _boundary_currents.end());

  // erase dirichlet nodes data structure
  _dirichlet_nodes.erase(_dirichlet_nodes.begin(),
      _dirichlet_nodes.end());

  // clear result vector
  _solution.erase(_solution.begin(), _solution.end());

  // clear variables vector
  _variables.erase(_variables.begin(), _variables.end());

  _rebuild_eq_system = true;
}

void
DriftDiffusion::solve(bool restart)
{
  prepare_solver();

  switch (_options.solver_method)
  {
    case GUMMEL: // TODO not implemented yet
      //solve_gummel(restart);
      //calculate_currents();
      break;
    default:
      try
      {
        solve_newton(restart);
        calculate_currents();
      }
      catch (PetscRuntimeError& e)
      {
        cerr << "Newton solver failed: " << e.get_reason() << "\n";
        throw (e);
      }
      break;
  }

  // aliases for nicer code
  EquationSystems& equation_systems = get_equation_system();
  
  // should perhaps only be done when the solution is requested
  equation_systems.build_solution_vector(_solution);
  for (int i = 0; i < _solution.size(); i++)
  {
    _solution[i] *= get_scaling().get_potential_scaling();
  }

}

double
DriftDiffusion::calculate_new_simulation_voltages(void)
{
  double step;
  double max_diff = 0;

  DD::Device::BoundaryList& boundaries = _device->get_boundaries();
  DD::Device::const_boundary_iterator it = boundaries.begin();
  const DD::Device::const_boundary_iterator end = boundaries.end();

  DD::Device::const_boundary_iterator max_id = end;

  for ( ; it != end; ++it)
  {
    BoundaryDescriptor* bd = *it;

    double diff = _simulation_voltages[bd] - _old_sim_voltages[bd];
    double abs_diff = fabs(diff);
    if (abs_diff > max_diff)
    {
      max_diff = abs_diff;
      step = diff;
      max_id = it;
    }
  }

  _simulation_voltages[*max_id] = _old_sim_voltages[*max_id] + step / 2;
  
  // TODO this is for test only
  cerr << "Try new voltage:\n";
  for ( it = boundaries.begin(); it != end; ++it)
  {
    BoundaryDescriptor* bd = *it;
    cerr << bd->get_id() << " : " << (_simulation_voltages[bd]) << " ";
  }
  cerr << "\n";


  return max_diff / 2.0;
}


void
DriftDiffusion::guess_equilibrium(NonlinearImplicitSystem& poisson) const
{
  
  const unsigned int u_var = poisson.variable_number("potential");
  const DofMap& dof_map_u = poisson.get_dof_map();
  vector<unsigned int> dof_indices_u;
  
  NumericVector<Number>& solution_u = *(poisson.solution);

  MeshBase::const_element_iterator el =
                                  get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  get_mesh().active_elements_end();

  const double phi0 = _scaling.get_potential_scaling();
  
  const unsigned int nn  = get_mesh().n_nodes();
  vector<unsigned short int> node_conn(nn);
  // Get the number of elements that share each node.  We will
  // compute the average value at each node.
  {
    vector<unsigned short int> node_conn_local(node_conn.size());
    
    MeshBase::const_element_iterator it =
      get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end =
      get_mesh().active_local_elements_end(); 

    for ( ; it != end; ++it)
      for (unsigned int n = 0; n < (*it)->n_nodes(); n++)
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

  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();
    SemiconductorModel* sc =
      _device->get_element_data().get_data(top_parent);
    dof_map_u.dof_indices(elem, dof_indices_u, u_var);
    for (int i = 0; i < elem->n_nodes(); i++)
    {
      solution_u.add(dof_indices_u[i],
          sc->get_equilibrium_fermi_level()
          / (phi0 * static_cast<Real>(node_conn[elem->node(i)])));
    }
  }
}

void
DriftDiffusion::set_solver_params(NonlinearSolver<Number>& solver,
    CalculationType calc_type)
{
  SolverClass& solver_class =
    static_cast<SolverClass&>(solver);
  
  SolverParameters& solver_params = _options.solver_params;

  const double phi0 = _scaling.get_potential_scaling();
  
  unsigned int nonlin_max_its = solver_params.nonlinear_max_iterations;
  // NOTE: we set nonlinear_max_iterations higher for equilibrium case
  //       because we assume that we need more to get the equilibrium
  //       solution than to get to a non-eq. solution
  if (calc_type == EQUILIBRIUM)
    nonlin_max_its *= 2;

  double sqrt_n = sqrt(get_mesh().n_nodes());
      
  solver_class.set_snes_options(solver_params.nonlinear_tolerance,
      solver_params.nonlinear_abs_tolerance * sqrt_n, nonlin_max_its,
      solver_params.ls_maxstep * sqrt_n / phi0);

  solver_class.set_ksp_options(solver_params.linear_tolerance,
      solver_params.linear_abs_tolerance * sqrt_n,
      solver_params.linear_max_iterations);

  solver_class.set_ksp_type(solver_params.ksp_type);
  solver_class.set_pc_type(solver_params.pc_type);
}


void
DriftDiffusion::solve_newton(bool restart)
{
  // aliases for nicer code
  Options& params = get_options();
  SolverParameters& solver_params = params.solver_params;
  
  DD::Device& device = *_device;
  Mesh& mesh = _device->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  EquationSystems& equation_systems = get_equation_system();
  Order approx_order = params.approximation_order;

  if (dim == 1)
    solver_params.ksp_type = KSPBCGS;
  else
    solver_params.ksp_type = KSPBCGSL;


  if (_rebuild_eq_system)
  {
    // the coupled DD system
    equation_systems.add_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

    NonlinearImplicitSystem& system =
      equation_systems.get_system<NonlinearImplicitSystem>(
          "drift-diffusion coupled");

    // we use PETSc
    system.nonlinear_solver =
      AutoPtr<NonlinearSolver<Number> >(new SolverClass);

    // set the options for the PETSc nonlinear solver
    set_solver_params(*system.nonlinear_solver, EQUILIBRIUM);

    system.add_variable("potential", approx_order);
    system.add_variable("fermi_e", approx_order);
    system.add_variable("fermi_h", approx_order);

    // we remember the equilibrium solution for future use
    system.add_vector("equilibrium solution");
    // we can remember a solution for future use
    system.add_vector("remembered solution");
    // for adaptive mesh refinement we need the old solution
    // befor a refinement step
    system.add_vector("old solution");

    system.nonlinear_solver->matvec = assemble<POISSON>;

    initialize_eq_system(equation_systems);

    // the simulation voltages for equilibrium are 0.0 ...
    _old_sim_voltages = _simulation_voltages;
    DD::Device::BoundaryList& boundaries = _device->get_boundaries();
    DD::Device::const_boundary_iterator it = boundaries.begin();
    for ( ; it != boundaries.end(); ++it)
    {
      _simulation_voltages[*it] = 0.0;
    }


    // make a rough guess of the equilibrium potential
    // TODO: the current guess (equil. potential) is not always a good thing
    //       so it's commented out
    guess_equilibrium(system);

    // try to solve the equilibrium potential
    //
    // If the SNES or KSP solvers diverge for some reason, we assume
    // that something is really bad and return. Only for the case that
    // it reached the maximum nonlinear iterations we give a second chance.
    try
    {
      system.solve();
      _n_nonlinear_iterations = system.n_nonlinear_iterations();
      _final_residual = system.final_nonlinear_residual();
    }
    catch (KSPDivergedError& e)
    {
      cerr << "*** PANIC: KSP diverged for equilibrium:\n";
      cerr << "   " << e.get_reason() <<
        " at iteration " << e.get_iteration() <<
        " (fnorm = " << e.get_fnorm() << ")\n";
      throw(e);
    }
    catch (SNESDivergedError& e)
    {
      cerr << "SNES diverged: " << e.get_reason() << "\n";
      if (e.get_reason() == SNES_DIVERGED_MAX_IT)
      {
        // we give him another chance to converge
        try
        {
          system.solve();
          _n_nonlinear_iterations = system.n_nonlinear_iterations();
          _final_residual = system.final_nonlinear_residual();
        }
        catch (PetscDivergedError& e)
        {
          cerr << "*** PANIC: SNES diverged for equilibrium:\n";
          cerr << "   " << e.get_reason() <<
            " at iteration " << e.get_iteration() <<
            " (fnorm = " << e.get_fnorm() << ")\n";
          throw(e);
        }
      }
      else
      {
        cerr << "*** PANIC: SNES diverged for equilibrium:\n";
        cerr << "   " << e.get_reason() <<
          " at iteration " << e.get_iteration() <<
          " (fnorm = " << e.get_fnorm() << ")\n";
        throw(e);
      }
    }

    system.get_vector("equilibrium solution") = *(system.solution);
    system.get_vector("equilibrium solution").close();

    system.get_vector("old solution") = *(system.solution);

    equation_systems.build_variable_names(_variables);

    // reset the simulation voltages
    _simulation_voltages = _old_sim_voltages;
    for (it = boundaries.begin(); it != boundaries.end(); ++it)
    {
      _old_sim_voltages[*it] = 0.0;
    }


    _rebuild_eq_system = false;
  }
  else
  {
    // restart from equilibrium solution
    if (restart)
    {
      NonlinearImplicitSystem& system =
        equation_systems.get_system<NonlinearImplicitSystem>(
            "drift-diffusion coupled");
    
      *(system.solution) = system.get_vector("equilibrium solution");

      // the simulation voltages for equilibrium are 0.0 ...
      DD::Device::BoundaryList& boundaries = _device->get_boundaries();
      DD::Device::const_boundary_iterator it = boundaries.begin();
      for ( ; it != boundaries.end(); ++it)
      {
        _old_sim_voltages[*it] = 0.0;
      }

      // just assign the old equilibrium solution doesn't work
      // very well. Why??
      if (device.get_mesh().n_elem() != device.get_mesh().n_active_elem())
        system.solve();
      //_simulation_voltages = _old_sim_voltages;
    }
  }

  NonlinearImplicitSystem& system =
    equation_systems.get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

  NumericVector<Number>& solution = *(system.solution);
  NumericVector<Number>& old_solution =
    system.get_vector("old solution");

  // set the solver parameters (they could have change since we made
  // the first calculation)
  set_solver_params(*system.nonlinear_solver, NONEQUILIBRIUM);


  // we don't simulate anything if the voltages didn't change
  if (_simulation_voltages == _old_sim_voltages) return;

  _options.C0_e = _options.n_max;
  _options.C0_h = _options.p_max;

  //system.nonlinear_solver->matvec = assemble<POISSON | ECURRENT>;
  system.nonlinear_solver->matvec = assemble<FULLYCOUPLED>;

  //
  // solve for the desired contact voltages
  //
  ContactData voltages = _simulation_voltages;

  bool reached = false;
  do
  {
    _simulation_voltages = voltages;
    bool retry = true;
    do
    {
      // try to solve the full step
      try
      {
        system.solve();
        
        retry = false;
      }
      // id it did not converge, so try half of the step
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
          cerr << " (Zero pivot during ILU.)";
        cerr << "\n";
      }

      if (retry)
      {
        double step = calculate_new_simulation_voltages();
        if (step < params.min_voltage_step) throw 1;

        // we have to clear the solver context, because the exception
        // left it in a unknown state such that future solves would fail.
        system.nonlinear_solver->clear();
        solution = old_solution;
      }
    } while (retry);

    _n_nonlinear_iterations = system.n_nonlinear_iterations();
    _final_residual = system.final_nonlinear_residual();
    _old_sim_voltages = _simulation_voltages;
    old_solution = solution;

    _options.C0_e = _options.n_max;
    _options.C0_h = _options.p_max;

    if (_simulation_voltages == voltages) reached = true;
  } while (!reached);


  if (params.mesh_refinement)
  {
    MeshRefinement mesh_refinement(device.get_mesh());
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
      
      find_dirichlet_nodes();

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


void
DriftDiffusion::solve_gummel(bool restart)
{
}


// TODO: currently not used
/*
void
DriftDiffusion::set_dirichlet_values(void)
{
  assert(_static_parameters.simulation_voltages->size() != 0);

  NonlinearImplicitSystem* poisson;
  NonlinearImplicitSystem* ecurrent;
  NonlinearImplicitSystem* hcurrent;

  switch (_options.solver_method)
  {
    case GUMMEL:
      poisson =
        &_eq_system->get_system<NonlinearImplicitSystem>("poisson");
      ecurrent =
        &_eq_system->get_system<NonlinearImplicitSystem>("ecurrent");
      hcurrent =
        &_eq_system->get_system<NonlinearImplicitSystem>("hcurrent");
      break;
    default:
      poisson = &_eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
      ecurrent = poisson;
      hcurrent = poisson;
      break;
  }

  ContactData& sim_voltages =
    *_static_parameters.simulation_voltages;

  // aliases for nicer code
  const Mesh& mesh = poisson->get_mesh();
  NumericVector<Number>& solution_u = *(poisson->solution);
  NumericVector<Number>& solution_en = *(ecurrent->solution);
  NumericVector<Number>& solution_ep = *(hcurrent->solution);

  const DofMap& dof_map_u = poisson->get_dof_map();
  const DofMap& dof_map_en = ecurrent->get_dof_map();
  const DofMap& dof_map_ep = hcurrent->get_dof_map();
  
  // numeric ids corresponding to the variables
  const unsigned int u_var = poisson->variable_number("potential");
  const unsigned int en_var = ecurrent->variable_number("fermi_e");
  const unsigned int ep_var = hcurrent->variable_number("fermi_h");

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // find all Dirichlet type boundaries (for at least one variable)
  set<BoundaryDescriptor*> dirichlet_boundaries;
  Device::BoundaryList& boundaries = _device->get_boundaries();
  Device::const_boundary_iterator it = boundaries.begin();
  const Device::const_boundary_iterator end = boundaries.end();
  for ( ; it != end; ++it)
  {
    if (((*it)->get_type("potential") == BoundaryDescriptor::DIRICHLET) ||
        ((*it)->get_type("fermi_e") == BoundaryDescriptor::DIRICHLET) ||
        ((*it)->get_type("fermi_h") == BoundaryDescriptor::DIRICHLET))
      dirichlet_boundaries.insert(*it);
  }

  set<BoundaryDescriptor*>::iterator not_dirichlet =
    dirichlet_boundaries.end();

  const BoundaryData& bound_data = _device->get_boundary_data();

  const BoundaryNodeList::const_iterator not_dirichlet_node =
          _dirichlet_nodes.end();
  BoundaryNodeList::const_iterator node_it;

  // loop over all level 0 boundary elements
  BoundaryData::const_iterator bd_it = bound_data.sides_begin();
  const BoundaryData::const_iterator bd_end = bound_data.sides_end();
  for ( ; bd_it != bd_end; ++bd_it)
  {
    if (dirichlet_boundaries.find(bd_it->second) != not_dirichlet)
    {

      // get the active family tree of this element
      vector<const Elem*> fam_tree;
      (bd_it->first).first->active_family_tree(fam_tree);

      vector<const Elem*>::const_iterator el = fam_tree.begin();
      for ( ; el != fam_tree.end(); ++el) 
      {
        const Elem* elem = *el;

        if (elem->on_boundary())
        {
          // get DOF indices
          dof_map_u.dof_indices(elem, dof_indices_u, u_var);
          dof_map_en.dof_indices(elem, dof_indices_en, en_var);
          dof_map_ep.dof_indices(elem, dof_indices_ep, ep_var);

          for (unsigned int i = 0; i < elem->n_nodes(); i++)
          {
            node_it = _dirichlet_nodes.find(elem->get_node(i));
            if (node_it != not_dirichlet_node)
            {
              const vector<double>* bc; 
              const BoundaryDescriptor* boundary_desc = node_it->second;
              if (boundary_desc->get_type("potential")
                  == BoundaryDescriptor::DIRICHLET)
              {
                bc = boundary_desc->get_coefficients("potential");
                solution_u.set(dof_indices_u[i],
                    (*bc)[2] / (*bc)[0] + sim_voltages[boundary_desc]);
              }

              if (boundary_desc->get_type("fermi_e")
                  == BoundaryDescriptor::DIRICHLET)
              {
                bc = node_it->second->get_coefficients("fermi_e");
                solution_en.set(dof_indices_en[i],
                    (*bc)[2] / (*bc)[0] - sim_voltages[boundary_desc]);
              }

              if (boundary_desc->get_type("fermi_h")
                  == BoundaryDescriptor::DIRICHLET)
              {
                bc = node_it->second->get_coefficients("fermi_h");
                solution_ep.set(dof_indices_ep[i],
                    (*bc)[2] / (*bc)[0] - sim_voltages[boundary_desc]);
              }
            }
          }
        }
      }
    }
  }
}
*/


double
DriftDiffusion::get_artificial_boundary_current(void)
{
  NonlinearImplicitSystem* poisson;
  NonlinearImplicitSystem* ecurrent;
  NonlinearImplicitSystem* hcurrent;

  switch (_options.solver_method)
  {
    case GUMMEL:
      poisson =
        &_eq_system->get_system<NonlinearImplicitSystem>("poisson");
      ecurrent =
        &_eq_system->get_system<NonlinearImplicitSystem>("ecurrent");
      hcurrent =
        &_eq_system->get_system<NonlinearImplicitSystem>("hcurrent");
      break;
    default:
      poisson = &_eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
      ecurrent = poisson;
      hcurrent = poisson;
      break;
  }


  // aliases for nicer code
  const Mesh& mesh = poisson->get_mesh();
  const NumericVector<Number>& solution_u = *(poisson->solution);
  const NumericVector<Number>& solution_en = *(ecurrent->solution);
  const NumericVector<Number>& solution_ep = *(hcurrent->solution);
  const DD::Device& device = *(_device);

  const DofMap& dof_map_u = poisson->get_dof_map();
  const DofMap& dof_map_en = ecurrent->get_dof_map();
  const DofMap& dof_map_ep = hcurrent->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  
  const double phi0 = _scaling.get_potential_scaling();
  
  double j0 = phi0;
  switch (dim)
  {
    case 1:
      j0 /= _options.mesh_units;
      break;
    case 3:
      j0 *= _options.mesh_units;
      break;
  }

  // numeric ids corresponding to the variables
  const unsigned int u_var = poisson->variable_number("potential");
  const unsigned int en_var = ecurrent->variable_number("fermi_e");
  const unsigned int ep_var = hcurrent->variable_number("fermi_h");
  
  FEType fe_type = poisson->variable_type(u_var);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qface(dim - 1, _options.integration_order);
  fe_face->attach_quadrature_rule(&qface);
  
  // Jacobian * quadrature weight at each integration point.   
  const vector<Real>& JxW = fe_face->get_JxW();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // this will hold the calculated properties of the semiconductor
  SemiconductorModel::CalculatedProperties properties;

  double tot_current = 0.0;
  
  MeshBase::const_element_iterator el =
                                  mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    // get DOF indices
    dof_map_u.dof_indices(elem, dof_indices_u, u_var);
    dof_map_en.dof_indices(elem, dof_indices_en, en_var);
    dof_map_ep.dof_indices(elem, dof_indices_ep, ep_var);

    SemiconductorModel* sc =
      device.get_element_data().get_data((*el)->top_parent());
    assert(sc != NULL);

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {
        int s_top = find_boundary(elem, s, top_parent);

        if (s_top != -1)
          continue; // is on a boundary

        fe_face->reinit(elem, s);

        int phi_size = phi.size();

        double current = 0.0;
        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          // get the solution value at the quadrature point
          Real u  = 0.0;
          Real en = 0.0;
          Real ep = 0.0;
          for (unsigned int i = 0; i < phi_size; i++)
          {
            u  += phi[i][qp] * solution_u(dof_indices_u[i]);
            en += phi[i][qp] * solution_en(dof_indices_en[i]);
            ep += phi[i][qp] * solution_ep(dof_indices_ep[i]);
          }

          // calculate densities etc.
          sc->calculate_all(u * phi0, en * phi0, ep * phi0, properties);
        
          Real n = properties.electron_density;
          Real p = properties.hole_density;

          Real mu_e = properties.electron_mobility;
          Real mu_h = properties.hole_mobility;

          Real cond_e = mu_e * n * Constants::e;
          Real cond_h = mu_h * p * Constants::e;

          for (unsigned int i = 0; i < phi_size; i++)
          {
            current += JxW[qp] * (dphi[i][qp] * face_normals[qp]) *
              (cond_e * solution_en(dof_indices_en[i]) +
               cond_h * solution_ep(dof_indices_ep[i]));
          }
        } // end loop over quadrature points

        tot_current += current;
      }
    } // end loop over elem sides
  } // end loop over elements

  return tot_current * j0;
}

void
DriftDiffusion::calculate_currents(void)
{

  // reset currents
  ContactData::iterator it =
    _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second = 0.0;

  NonlinearImplicitSystem* system =
    &_eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");


  // aliases for nicer code
  const Mesh& mesh = system->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);
  const DD::Device& device = *(_device);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = _scaling.get_potential_scaling();

  // the scaling for the current
  // NOTE: e current eg. is integral(Jacobian * mu * n * nabla Ef * face_normal)
  //       - Ef has to be scaled back by phi0
  //       - nabla has to be scaled back by 1 / x0, but the mesh is not scaled
  //         (apart from the factor 'mesh_units')
  double j0 = phi0;
  switch (dim)
  {
    case 1:
      j0 /= _options.mesh_units;
      break;
    case 3:
      j0 *= _options.mesh_units;
      break;
  }
  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = _options.integration_order;
  
  QGauss qface(dim - 1, integration_order);
  fe_face->attach_quadrature_rule(&qface);

  
  
  // Jacobian * quadrature weight at each integration point.   
  const vector<Real>& JxW = fe_face->get_JxW();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // this will hold the calculated properties of the semiconductor
  SemiconductorModel::CalculatedProperties properties;

  MeshBase::const_element_iterator el =
                                  mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    SemiconductorModel* sc =
      device.get_element_data().get_data((*el)->top_parent());
    assert(sc != NULL);

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {
        int s_top = find_boundary(elem, s, top_parent);

        if (s_top == -1)
          continue; // is not on a boundary

        BoundaryData::ElementSide side(top_parent, s_top);
        BoundaryDescriptor* boundary_desc =
          device.get_boundary_data().get_data(side);


        // only for dim > 1 we need to integrate
        if (dim > 1)
        {
          fe_face->reinit(elem, s);

          int phi_size = phi.size();

          double current = 0.0;

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // get the solution value at the quadrature point
            Real u  = 0.0;
            Real en = 0.0;
            Real ep = 0.0;
            Real dEfn = 0.0;
            Real dEfp = 0.0;
            for (unsigned int i = 0; i < phi_size; i++)
            {
              u  += phi[i][qp] * solution(dof_indices_u[i]);
              en += phi[i][qp] * solution(dof_indices_en[i]);
              ep += phi[i][qp] * solution(dof_indices_ep[i]);
              double tmp = dphi[i][qp] * face_normals[qp];
              dEfn += tmp * solution(dof_indices_en[i]);
              dEfp += tmp * solution(dof_indices_ep[i]);
            }

            // calculate densities etc.
            sc->calculate_all(u * phi0, en * phi0, ep * phi0, properties);

            Real n = properties.electron_density;
            Real p = properties.hole_density;

            Real mu_e = properties.electron_mobility;
            Real mu_h = properties.hole_mobility;

            Real cond_e = mu_e * n * Constants::e;
            Real cond_h = mu_h * p * Constants::e;

            current += JxW[qp] * (cond_e * dEfn + cond_h * dEfp);
          } // end loop over quadrature points

          _boundary_currents[boundary_desc] += current;
        }
        else
        {
          vector<Point> v(1, elem->point(s));
          fe_face->reinit(elem, &v);

          double current = 0.0;

          // s is the node of the element lying on the boundary
          Real u  = solution(dof_indices_u[s]);
          Real en = solution(dof_indices_en[s]);
          Real ep = solution(dof_indices_ep[s]);

          Real dEfn = 0.0;
          Real dEfp = 0.0;
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            dEfn += dphi[n][0](0) * solution(dof_indices_en[n]);
            dEfp += dphi[n][0](0) * solution(dof_indices_ep[n]);
          }

          // calculate densities etc.
          sc->calculate_all(u * phi0, en * phi0, ep * phi0, properties);

          Real n = properties.electron_density;
          Real p = properties.hole_density;

          Real mu_e = properties.electron_mobility;
          Real mu_h = properties.hole_mobility;

          Real cond_e = mu_e * n * Constants::e;
          Real cond_h = mu_h * p * Constants::e;

          _boundary_currents[boundary_desc] =
            (cond_e * dEfn + cond_h * dEfp);
        }
      }
    } // end loop over elem sides
  } // end loop over elements

  // scale the current to normal units
  it = _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second *= j0;
}


// implementation taken from libmesh equation_systems.C
  void
DriftDiffusion::build_densities(vector<double>& densities,
    vector<string>& names)
{
  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_nodes();

  const unsigned int n_vars  = 3;
  names.resize(n_vars);
  names[0] = "e_density";
  names[1] = "h_density";
  names[2] = "recomb_rate";

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

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // this will hold the calculated properties of the semiconductor
  SemiconductorModel::CalculatedProperties properties;

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
      dof_map.dof_indices(elem, dof_indices_en, en_var);
      dof_map.dof_indices(elem, dof_indices_ep, ep_var);

      SemiconductorModel* sc =
        device.get_element_data().get_data(elem->top_parent());
      assert(sc != NULL);

      assert(elem->n_nodes() == dof_indices_u.size());

      //sol.resize(dof_indices.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        
        sc->calculate_all(u, en, ep, properties);

        assert (node_conn[elem->node(n)] != 0);

        unsigned int id = n_vars * elem->node(n);
        double nodal_val = properties.electron_density;
        local[id] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = properties.hole_density;
        local[id + 1] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = properties.net_electron_recombination_rate;
        local[id + 2] +=
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


//template <DriftDiffusion::EquationType T, int dim>
template <int T>
void
DriftDiffusion::assemble(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian)
{
  
  // references for nicer code
  const Mesh& mesh = _this->get_mesh();
  EquationSystems& eq_sys = _this->get_equation_system();
  NonlinearImplicitSystem& system =
    eq_sys.get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

  const unsigned int dim = mesh.mesh_dimension();
  
  const DD::Device& device = _this->get_device();
  const Options& params = _this->get_options();
  Options& options = _this->get_options();

  ContactData& simulation_voltages = _this->_simulation_voltages;
  BoundaryNodeList& dirichlet_nodes = _this->_dirichlet_nodes;


  //
  // some scaling stuff...
  // 
  // NOTE: the mesh and all paramters were not explicitly scaled, so
  //       we have to treat scaling by explicit division/multiplication
  //       
  // density scaling for electrons
  double C0_e = options.C0_e;
  // density scaling for holes
  double C0_h = options.C0_h;
  // maximum density of electrons
  double n_max = options.n_max;
  // maximum density of holes
  double p_max = options.p_max;
  // the scaling parameters
  const Scaling& scaling = _this->get_scaling();
  // the scaling parameter for the poisson eq.
  // The factor 1e-2 comes from the fact, that we are calculating in cm!
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  // x 1e4 because we calculate in cm
  const double P0 = (Constants::e * x0 * C0) * 1e4;
  const double R0_e = C0_e / scaling.get_time_scaling();
  const double R0_h = C0_h / scaling.get_time_scaling();
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
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_e");
  const unsigned int ep_var = system.variable_number("fermi_h");
  
  FEType fe_type = system.variable_type(u_var);

  // the finite element
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, params.integration_order);
  fe->attach_quadrature_rule(&qrule);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = params.integration_order;
  
  QGauss qface(dim - 1, integration_order);
  fe_face->attach_quadrature_rule(&qface);

  
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
    Kuu(Ke), Kun(Ke), Kup(Ke),
    Knu(Ke), Knn(Ke), Knp(Ke),
    Kpu(Ke), Kpn(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fn(Fe),
    Fp(Fe);

  DenseSubVector<Number>
    Xu(X),
    Xn(X),
    Xp(X);


  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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
    Elem* elem_tmp = *el;
    const Elem* top_parent = (*el)->top_parent();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    unsigned int n_dofs     = dof_indices_u.size();
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
    //        | Kuu Kun Kup |        | Fu |
    //   Ke = | Knu Knn Knp |;  Fe = | Fn |
    //        | Kpu Kpn Kpp |        | Fp |
    //         -           -          -  -
    //
    Kuu.reposition(0, 0, n_dofs, n_dofs);
    Kun.reposition(0, n_dofs, n_dofs, n_dofs);
    Kup.reposition(0, 2 * n_dofs, n_dofs, n_dofs);
    //
    Knu.reposition(n_dofs, 0, n_dofs, n_dofs);
    Knn.reposition(n_dofs, n_dofs, n_dofs, n_dofs);
    Knp.reposition(n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    //
    Kpu.reposition(2 * n_dofs, 0, n_dofs, n_dofs);
    Kpn.reposition(2 * n_dofs, n_dofs, n_dofs, n_dofs);
    Kpp.reposition(2 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    //
    Fu.reposition(0, n_dofs);
    Fn.reposition(n_dofs, n_dofs);
    Fp.reposition(2 * n_dofs, n_dofs);
    //
    Xu.reposition(0, n_dofs);
    Xn.reposition(n_dofs, n_dofs);
    Xp.reposition(2 * n_dofs, n_dofs);


    DriftDiffusionProperties* sc =
      device.get_element_data().get_data(top_parent);
    assert(sc != NULL);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the solution values at the quadrature points
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * Xu(i);
        en += phi[i][qp] * Xn(i);
        ep += phi[i][qp] * Xp(i);
      }

      // calculate densities etc.
      sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, q_point[qp], elem);
      double n  = sc->get_electron_density();
      double p  = sc->get_hole_density();
      double epsilon = sc->get_relative_permittivity();
      double l2_eps = l2 * epsilon;
      
      // remember the maximum densities
      n_max = (n_max > n) ? n_max : n;
      p_max = (p_max > p) ? p_max : p;

      // NOTE: sigma_e = mu_e * n is the electron conductivity
      Real sigma_e = sc->get_electron_conductivity() / (mu0 * C0_e);
      Real sigma_h = sc->get_hole_conductivity() / (mu0 * C0_h);
    
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
          Real laplace =
            -J * (dphi[i][qp] * dphi[j][qp]) * x0_mesh * x0_mesh;
          
          if (T & POISSON)
            Kuu(i,j) += l2_eps * laplace;
          
          if (T & ECURRENT)
              Knn(i,j) += sigma_e * laplace;
          
          if (T & HCURRENT)
              Kpp(i,j) += sigma_h * laplace;
        }

        if (!(T & POISSON))
          Kuu(i,i) += 1;
        
        if (!(T & ECURRENT))
          Knn(i,i) += 1;
        
        if (!(T & HCURRENT))
          Kpp(i,i) += 1;
      }
            
      // 
      // for jacobian compute the other contributions
      // 
      if (jacobian != NULL)
      {
        double drho[3];
        double dRn[3];
        double dRp[3];
        for (int id = 0; id < 3; id++)
        {
          drho[id] = phi0 / C0 * sc->get_charge_density_derivatives()[id];
          dRn[id] = phi0 / R0_e
            * sc->get_net_electron_recombination_rate_derivatives()[id];
          dRp[id] = phi0 / R0_h
            * sc->get_net_hole_recombination_rate_derivatives()[id];
        }

        // d(sigma_n)/du * element-jacobian
        // sigma_n = mu_n * n means the conductivity of electrons
        Real dsigma_e = J * phi0 / (mu0 * C0_e)
          * sc->get_electron_conductivity_derivatives()[0];
        Real dsigma_h = J * phi0 / (mu0 * C0_h)
          * sc->get_hole_conductivity_derivatives()[0];

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            
            Real dsigma_e_x_phi = dsigma_e * phi[j][qp];
            Real dsigma_h_x_phi = dsigma_h * phi[j][qp];
            for (unsigned int k = 0; k < n_dofs; k++)
            {
              Real laplace = -(dphi[i][qp] * dphi[k][qp]) * x0_mesh * x0_mesh;
            
              if (T & ECURRENT)
              {
                Real elem_contrib =
                  dsigma_e_x_phi * laplace * Xn(k);

                if (T & POISSON)
                  Knu(i,j) += elem_contrib;

                Knn(i,j) += elem_contrib;
              }

              if (T & HCURRENT)
              {
                Real elem_contrib =
                  dsigma_h_x_phi * laplace * Xp(k);
                
                if (T & POISSON)
                  Kpu(i,j) += elem_contrib;

                Kpp(i,j) += elem_contrib;
              }
            }

            // The dFe_i/dX_j part
            Real phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (T & POISSON)
            {
              Kuu(i,j) += drho[0] * phi_i_x_phi_j;
              
              if (T & ECURRENT)
                Kun(i,j) += drho[1] * phi_i_x_phi_j;

              if (T & HCURRENT)
                Kup(i,j) += drho[2] * phi_i_x_phi_j;
            }            
            
            if (T & ECURRENT)
            {
              if (T & POISSON)
                Knu(i,j) -= dRn[0] * phi_i_x_phi_j;

              Knn(i,j) -= dRn[1] * phi_i_x_phi_j;

              if (T & HCURRENT)
                Knp(i,j) -= dRn[2] * phi_i_x_phi_j;
            }

            if (T & HCURRENT)
            {
              if (T & POISSON)
                Kpu(i,j) += dRp[0] * phi_i_x_phi_j;

              if (T & ECURRENT)
                Kpn(i,j) += dRp[1] * phi_i_x_phi_j;
              
              Kpp(i,j) += dRp[2] * phi_i_x_phi_j;
            }

          }
        }
      }

      // if we are doing residual, calculate rhs contribution (i.e. -Fe)
      if (residual != NULL)
      {
        // charge density
        Real J_x_rho = J * sc->get_charge_density() / C0;
        Real J_x_P0 = J / P0;

        // net recombination rate
        Real J_x_Rn = J * sc->get_net_electron_recombination_rate() / R0_e;
        Real J_x_Rp = J * sc->get_net_hole_recombination_rate() / R0_h;

        RealVectorValue P(sc->get_total_polarization());
        P *= J_x_P0;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          Real net_recomb_e = J_x_Rn * phi[i][qp];
          Real net_recomb_h = J_x_Rp * phi[i][qp];
          
          if (T & POISSON)
            Fu(i) += J_x_rho * phi[i][qp] + (P * dphi[i][qp]) * x0_mesh;
          else
            Fu(i) -= Xu(i);
          
          if (T & ECURRENT)
            Fn(i) -= net_recomb_e;
          else
            Fn(i) -= Xn(i);

          if (T & HCURRENT)
            Fp(i) += net_recomb_h;
          else
            Fp(i) -= Xp(i);
        }
      }

    } // end loop over quadrature points

    // now loop over the element sides to find boundary elements
    // and to include von Neumann and mixed type boundary conditions
    // 
    // NOTE 1:
    // we dont apply BC for nabla(Ef) but for the particle
    // flux mu * n * nabla(Ef)
    //
    // NOTE 2:
    // 1D case needs special treatment as 0D boundary elements do not
    // exist in libmesh...
    // 
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {
        int s_top = find_boundary(elem, s, top_parent);

        BoundaryData::ElementSide side(top_parent, s_top);
        BoundaryDescriptor* boundary_desc =
          device.get_boundary_data().get_data(side);

        vector<double> coeff(3, 0.0);
        vector<double> value(3, 0.0);
        if (boundary_desc != NULL)
        {
          const vector<double>* coeffs;
          double c, v;

          if (T & POISSON)
          {
            if (boundary_desc->get_type("potential") !=
                BoundaryDescriptor::DIRICHLET)
            {
              coeffs =  boundary_desc->get_coefficients("potential");
              assign_boundary_values(c, v, *coeffs);
              coeff[0] = c * x0;
              value[0] = v * x0 / phi0;
            }
          }
          if (T & ECURRENT)
          {
            if (boundary_desc->get_type("fermi_e") !=
                BoundaryDescriptor::DIRICHLET)
            {
              coeffs = boundary_desc->get_coefficients("fermi_e");
              assign_boundary_values(c, v, *coeffs);
              coeff[1] = c * x0;
              value[1] = v * x0 / phi0;
            }
          }
          if (T & HCURRENT)
          {
            if (boundary_desc->get_type("fermi_h") !=
                BoundaryDescriptor::DIRICHLET)
            {
              coeffs = boundary_desc->get_coefficients("fermi_h");
              assign_boundary_values(c, v, *coeffs);
              coeff[2] = c * x0;
              value[2] = v * x0 / phi0;
            }
          }
        }

        //
        // NOTE: we have to integrate over the boundary also if there are
        //       no contacts because there could be polarization.
        //
        const vector<vector<Real> >&  phi_face =
          fe_face->get_phi();
          
        // physical coordinates of the quadrature points
        const vector<Point>& q_point = fe_face->get_xyz();

        const vector<Point>& face_normals = fe_face->get_normals();

        const vector<Real>& JxW_face = fe_face->get_JxW();

        if (dim > 1)
        {
          fe_face->reinit(elem, s);

          int phi_size = phi_face.size();

          // now integrate to include von Neumann and mixed type BCs
          // and polarization
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {

            // get the solution values at the quadrature points
            Real u  = 0.0;
            Real en = 0.0;
            Real ep = 0.0;
            for (unsigned int i = 0; i < n_dofs; i++)
            {
              u  += phi_face[i][qp] * Xu(i);
              en += phi_face[i][qp] * Xn(i);
              ep += phi_face[i][qp] * Xp(i);
            }

            // calculate densities etc.
            sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, q_point[qp], elem);
            double epsilon = sc->get_relative_permittivity();
            double l2_eps = l2 * epsilon;

            // the jacobian x weight x scaling
            double J = JxW_face[qp] / Jface_scale;

            // first the contributions to Ke_ij
            for (unsigned int i = 0; i < n_dofs; i++)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {

                Real phi_i_x_phi_j =
                  J * phi_face[i][qp] * phi_face[j][qp];

                if (T & POISSON)
                  Kuu(i,j) -= l2_eps * coeff[0] * phi_i_x_phi_j;

                if (T & ECURRENT)
                  Knn(i,j) -= coeff[1] * phi_i_x_phi_j;

                if (T & HCURRENT)
                  Kpp(i,j) -= coeff[2] * phi_i_x_phi_j;
              }
            }

            // contribution to -Fe_i
            if (residual != NULL)
            {
              RealVectorValue P(sc->get_total_polarization());
              double Pn = (P * face_normals[qp]) / P0;
              double value_u = J * (l2_eps * value[0] - Pn);
              double value_n = J * value[1];
              double value_p = J * value[2];

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                if (T & POISSON)
                  Fu(i) += value_u * phi_face[i][qp];

                if (T & ECURRENT)
                  Fn(i) += value_n * phi_face[i][qp];

                if (T & HCURRENT)
                  Fp(i) += value_p * phi_face[i][qp];
              }
            } 
          }
        }
        else // i.e. dim == 1
        {
          // s is the node of the element lying on the boundary
          Real u  = Xu(s);
          Real en = Xn(s);
          Real ep = Xp(s);

          // calculate densities etc.
          sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep,
            elem->point(s), elem);
          double epsilon = sc->get_relative_permittivity();
          double l2_eps = l2 * epsilon;

          // first the contributions to Ke_ij
          if (T & POISSON)
            Kuu(s,s) -= l2_eps * coeff[0];

          if (T & ECURRENT)
            Knn(s,s) -= coeff[1];

          if (T & HCURRENT)
            Kpp(s,s) -= coeff[2];

          // contribution to -Fe_i
          if (residual != NULL)
          {
            double Pn =  sc->get_total_polarization()(0) / P0;
            double value_u = l2_eps * value[0] - Pn;
            double value_n = value[1];
            double value_p = value[2];

            if (T & POISSON)
              Fu(s) += value_u;

            if (T & ECURRENT)
              Fn(s) += value_n;

            if (T & HCURRENT)
              Fp(s) += value_p;
          }
        }
      }
    } // end loop over element sides



    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    // NOTE: this changes dof_indices that's why the application of
    //       Dirichlet type BCs needs special care
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);


    //
    // now as last thing we apply Dirichlet type Bcs
    //
    //double l2_eps0 = l2 * Constants::e0;
    BoundaryNodeList::const_iterator node_it;
    const BoundaryNodeList::const_iterator end =
      dirichlet_nodes.end();
    if (Ke.m() == n_dofs_tot)
    {
      // no constrained nodes, so everything is easy

      // loop over all nodes and check if it is a dirichlet type node
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        node_it = dirichlet_nodes.find(elem->get_node(i));
        if (node_it != end)
        {
          const BoundaryDescriptor* boundary_desc = node_it->second;

          if (T & POISSON)
          {
            if (boundary_desc->get_type("potential") ==
                BoundaryDescriptor::DIRICHLET)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {
                Kuu(i,j) = 0.0; Kun(i,j) = 0.0; Kup(i,j) = 0.0;
              }

              const vector<double>* bc =
                boundary_desc->get_coefficients("potential");
              // TODO what is better???
              //Kuu(i,i) = -l2_eps0;
              //Fu(i)    = l2_eps0 * ((*bc)[2] / (*bc)[0] +
              Kuu(i,i) = 1;
              Fu(i)    = -((*bc)[2] / (*bc)[0] +
                  simulation_voltages[boundary_desc]) / phi0;
            }
          }

          if (T & ECURRENT)
          {
            if (boundary_desc->get_type("fermi_e") ==
                BoundaryDescriptor::DIRICHLET)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {
                Knu(i,j) = 0.0; Knn(i,j) = 0.0; Knp(i,j) = 0.0;
              }

              const vector<double>* bc =
                node_it->second->get_coefficients("fermi_e");
              Knn(i,i) = 1;
              Fn(i)    = -((*bc)[2] / (*bc)[0] -
                  simulation_voltages[boundary_desc]) / phi0;
            }
          }

          if (T & HCURRENT)
          {
            if (boundary_desc->get_type("fermi_h") ==
                BoundaryDescriptor::DIRICHLET)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {
                Kpu(i,j) = 0.0; Kpn(i,j) = 0.0; Kpp(i,j) = 0.0;
              }

              const vector<double>* bc =
                node_it->second->get_coefficients("fermi_h");
              Kpp(i,i) = 1;
              Fp(i)    = -((*bc)[2] / (*bc)[0] -
                  simulation_voltages[boundary_desc]) / phi0;
            }
          }
        }
      }
    }
    else
    {
      // Some nodes are constrained, so we have messed up our
      // matrix and vector. In particular, we could have included
      // nodes on Dirichlet boundaries.
      // We will look for them on the parent element(s) to apply
      // proper boundary conditions

      n_dofs_tot = dof_indices.size();

      // it's possible, that a node of the parent element is
      // also a hanging node. In this case we have to look at the
      // grand parent
      bool is_done = false;
      const Elem* parent;
      while (!is_done)
      {
        is_done = true;
        parent = elem->parent();
        elem = parent;

        assert(parent != NULL);

        dof_map.dof_indices(parent, dof_indices_u, u_var);
        dof_map.dof_indices(parent, dof_indices_en, en_var);
        dof_map.dof_indices(parent, dof_indices_ep, ep_var);

        // loop over the nodes of the parent element
        unsigned int n_nodes = parent->n_nodes();
        for (unsigned int i = 0; i < n_nodes; i++)
        {
          //if (dof_map.is_constrained_dof(dof_indices_u[i]))
          //  is_done = false;

          node_it = dirichlet_nodes.find(parent->get_node(i));
          if (node_it != end)
          {
            const BoundaryDescriptor* boundary_desc = node_it->second;

            // loop over all DOFs occurring in the constrained matrix
            for (unsigned int id = 0; id < n_dofs_tot; id++)
            {

              if (T & POISSON)
              {
                if (boundary_desc->get_type("potential") ==
                    BoundaryDescriptor::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_u[i])
                  {
                    for (unsigned int j = 0; j < n_dofs_tot; j++)
                      Ke(id,j) = 0.0;

                    const vector<double>* bc =
                      boundary_desc->get_coefficients("potential");
                    Ke(id,id) = 1;
                    Fe(id)    = -((*bc)[2] / (*bc)[0]
                        + simulation_voltages[boundary_desc]) / phi0;
                  }
                }
              }

              if (T & ECURRENT)
              {
                if (boundary_desc->get_type("fermi_e") ==
                    BoundaryDescriptor::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_en[i])
                  {
                    for (unsigned int j = 0; j < n_dofs_tot; j++)
                      Ke(id,j) = 0.0;

                    const vector<double>* bc =
                      boundary_desc->get_coefficients("fermi_e");
                    Ke(id,id) = 1;
                    Fe(id)    = -((*bc)[2] / (*bc)[0]
                        - simulation_voltages[boundary_desc]) / phi0;
                  }
                }
              }

              if (T & HCURRENT)
              {
                if (boundary_desc->get_type("fermi_h") ==
                    BoundaryDescriptor::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_ep[i])
                  {
                    for (unsigned int j = 0; j < n_dofs_tot; j++)
                      Ke(id,j) = 0.0;

                    const vector<double>* bc =
                      boundary_desc->get_coefficients("fermi_h");
                    Ke(id,id) = 1;
                    Fe(id)    = -((*bc)[2] / (*bc)[0]
                        - simulation_voltages[boundary_desc]) / phi0;
                  }
                }
              }

            } // end loop over all DOFs 
          }
        } // end loop over the nodes of the parent element
      }
    }

    if (residual != NULL)
    {
      // for residual, add K*X to -Fe, as R = Ke*X - Fe
      for (unsigned int i = 0; i < n_dofs_tot; i++)
        for (unsigned int j = 0; j < n_dofs_tot; j++)
          Fe(i) += Ke(i,j) * x(dof_indices[j]);

      residual->add_vector(Fe, dof_indices);
    }
    else
      jacobian->add_matrix(Ke, dof_indices);

  } // end loop over elements

  // put the maximum densities back into the options
  options.n_max = n_max;
  options.p_max = p_max;
} 

