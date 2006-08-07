// $Id$

// module includes
#include "DriftDiffusion.h"
#include "DDevice.h"
#include "Scaling.h"
#include "ElementData.h"
#include "BoundaryData.h"
#include "ElectricalContact.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"
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
using namespace DriftDiffusionDefs;


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
    solver_method(NEWTON),
    max_gummel_iterations(5),
    mesh_units(1e-4), // default mesh units are um
    scaling_type(Scaling::UNITS),
    coupling(POISSON),
    n_max(1), p_max(1),
    C0_e(1), C0_h(1),
    linearize_continuity_eq(false),
    artificial_drift(false),
    local_scaling(false)
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
    solver_method(rhs.solver_method),
    max_gummel_iterations(rhs.max_gummel_iterations),
    solver_params(rhs.solver_params),
    mesh_units(rhs.mesh_units),
    scaling_type(rhs.scaling_type),
    coupling(rhs.coupling),
    n_max(rhs.n_max),
    p_max(rhs.p_max),
    C0_e(rhs.C0_e),
    C0_h(rhs.C0_h),
    linearize_continuity_eq(rhs.linearize_continuity_eq),
    artificial_drift(rhs.artificial_drift),
    local_scaling(rhs.local_scaling)
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
    solver_method = rhs.solver_method;
    max_gummel_iterations = rhs.max_gummel_iterations;
    solver_params = rhs.solver_params;
    mesh_units = rhs.mesh_units;
    scaling_type = rhs.scaling_type;
    coupling = rhs.coupling;
    n_max = rhs.n_max;
    p_max = rhs.p_max;
    C0_e = rhs.C0_e;
    C0_h = rhs.C0_h;
    linearize_continuity_eq = rhs.linearize_continuity_eq;
    artificial_drift = rhs.artificial_drift;
    local_scaling = rhs.local_scaling;
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
    pc_type(PCILU)
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
  cleanup_solver();
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
  double phi0 = SimulationOptions::T * Constants::k_B;
  double mu0 = -1;
  double C0 = 1e-12;
  double ni0 = 1e-12; // let 1 be the minimum for ni0
  double eps0 = -1;
  
  // find minimum or maximum by looping over all elements
  MeshBase::const_element_iterator el =
                                  get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  get_mesh().active_local_elements_end();
  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    DriftDiffusionProperties* sc =
      _device->get_element_data().get_data(top_parent);

    sc->reinit(elem);
    sc->calculate_all(sc->get_equilibrium_fermi_level(), 0.0, 0.0,
        elem->centroid());
    
    double mu = sc->get_hole_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;
    mu = sc->get_electron_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;

    //double C = fabs(sc->get_net_doping_density());
    double C = fabs(sc->get_ionized_donor_density() -
        sc->get_ionized_acceptor_density());
    C0 = (C0 > C) ? C0 : C;

    double ni = sc->get_intrinsic_density();
    ni0 = (ni0 > ni) ? ni0 : ni;

    double eps = sc->get_relative_permittivity();
    eps0 = (eps0 > eps) ? eps0 : eps;
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
  //_scaling.set_potential_scaling(1.0);
  _scaling.set_length_scaling(x0 * _options.mesh_units);
  _scaling.set_mobility_scaling(mu0);
  _scaling.set_density_scaling(C0);
}


void
DriftDiffusion::set_simulation_voltage(const string& boundary,
    double voltage)
{
  ElectricalContact* desc = _device->get_boundary(boundary);

  if (desc != NULL)
    _simulation_voltages[desc] = voltage;
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

  system.get_vector("old solution") = *(system.solution);

}

void
DriftDiffusion::set_electron_fermi_level(double Ef_n)
{
  NonlinearImplicitSystem& system =
    _eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
  
  const unsigned int var = system.variable_number("fermi_e");
  const double phi0 = _scaling.get_potential_scaling();
  double level = -Ef_n / phi0;

  Mesh& mesh = get_mesh();
  Mesh::node_iterator it = mesh.active_nodes_begin();
  const Mesh::node_iterator end = mesh.active_nodes_end();

  for ( ; it != end; ++it)
  {
    const Node* node = *it;
    unsigned int id = node->dof_number(system.number(), var, 0);
    system.solution->set(id, level);
  }
}

void
DriftDiffusion::set_hole_fermi_level(double Ef_p)
{
  NonlinearImplicitSystem& system =
    _eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
  
  const unsigned int var = system.variable_number("fermi_h");
  const double phi0 = _scaling.get_potential_scaling();
  double level = -Ef_p / phi0;

  Mesh& mesh = get_mesh();
  Mesh::node_iterator it = mesh.active_nodes_begin();
  const Mesh::node_iterator end = mesh.active_nodes_end();

  for ( ; it != end; ++it)
  {
    const Node* node = *it;
    unsigned int id = node->dof_number(system.number(), var, 0);
    system.solution->set(id, level);
  }

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
  set<ElectricalContact*> dirichlet_boundaries;

  DD::Device::BoundaryList& boundaries = _device->get_boundaries();
  DD::Device::const_boundary_iterator it = boundaries.begin();
  const DD::Device::const_boundary_iterator end = boundaries.end();
  for ( ; it != end; ++it)
  {
    if (((*it)->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
        || ((*it)->get_type(POTENTIAL) == ElectricalContact::PINNING)
        || ((*it)->get_type(FERMIE) == ElectricalContact::DIRICHLET)
        || ((*it)->get_type(FERMIH) == ElectricalContact::DIRICHLET))
      dirichlet_boundaries.insert(*it);
  }

  set<ElectricalContact*>::iterator not_dirichlet =
    dirichlet_boundaries.end();

  // loop over all level 0 boundary elements
  BoundaryData::const_iterator bd_it = bound_data.sides_begin();
  const BoundaryData::const_iterator bd_end = bound_data.sides_end();
  for ( ; bd_it != bd_end; ++bd_it)
  {
    if (dirichlet_boundaries.find(bd_it->second) != not_dirichlet)
    {
      const Elem* elem = (bd_it->first).first;

      // dont't go on if it is an inactive element
      if (elem->refinement_flag() == Elem::INACTIVE) continue;

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
DriftDiffusion::reset_solver(void)
{
  if (!_rebuild_eq_system)
  {
    assert(_eq_system != NULL);
    _eq_system->delete_system("drift-diffusion coupled");
    _rebuild_eq_system = true;
  }
}

void
DriftDiffusion::cleanup_solver(void)
{
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

  reset_solver();
}

void
DriftDiffusion::solve(void)
{

  assert(_rebuild_eq_system == false);

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;

  switch (_options.solver_method)
  {
    case GUMMEL:
      try { solve_gummel(); }
      catch (PetscRuntimeError& e)
      {
        cerr << "Gummel iteration failed: " << e.get_reason() << "\n";
        throw (e);
      }
      calculate_currents();
      break;
    default: // Newton method
      try { solve_newton(); }
      catch (PetscRuntimeError& e)
      {
        cerr << "Newton solver failed: " << e.get_reason() << "\n";
        throw (e);
      }
      calculate_currents();
      break;
  }

  build_solution_vector(_solution);
  
  // update the list which contains all elements of this simulation
  update_element_list();
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
    ElectricalContact* bd = *it;

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
    ElectricalContact* bd = *it;
    cerr << bd->get_id() << " : " << (_simulation_voltages[bd]) << " ";
  }
  cerr << "\n";


  return max_diff / 2.0;
}


void
DriftDiffusion::guess_equilibrium(void)
{
  
  NonlinearImplicitSystem& poisson =
    get_equation_system().get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

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
    DriftDiffusionProperties* sc =
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
DriftDiffusion::update_element_list(void)
{

  MeshBase::const_element_iterator el =
                                  get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  get_mesh().active_elements_end();

  _element_list.clear();
  for ( ; el != end_el; ++el)
  {
    _element_list.insert(*el);
  }
}

void
DriftDiffusion::init(void)
{
  if (!_rebuild_eq_system) return;

  EquationSystems& equation_systems = get_equation_system();

  // the coupled DD system
  NonlinearImplicitSystem& system =
    equation_systems.add_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

  // we use PETSc
  system.nonlinear_solver =
    AutoPtr<NonlinearSolver<Number> >(new SolverClass);

  system.add_variable("potential", libMeshEnums::FIRST);
  system.add_variable("fermi_e", libMeshEnums::FIRST);
  system.add_variable("fermi_h", libMeshEnums::FIRST);

  // we can remember a solution for future use
  system.add_vector("remembered solution");
  // for adaptive mesh refinement we need the old solution
  // befor a refinement step
  system.add_vector("old solution");

  // for local scaling
  // TODO will probably be deleted
  system.add_vector("scaling");


  // set some parameters, which we don't use in this way, though.
  // But they are needed by libmesh
  SolverParameters& solver_params =
    get_options().solver_params;

  equation_systems.parameters.set<unsigned int>(
    "nonlinear solver maximum iterations") = 
      solver_params.nonlinear_max_iterations;
  
  equation_systems.parameters.set<Real>("nonlinear solver tolerance") =
    solver_params.nonlinear_tolerance;


  // finally initialize the newly created system
  system.init();

  // compute the scaling factors
  compute_scaling(get_options().scaling_type);
  
  // make a rough guess of the equilibrium potential
  guess_equilibrium();

  // prepare a list for the boundary voltages
  DD::Device::BoundaryList& boundaries = _device->get_boundaries();
  DD::Device::const_boundary_iterator it = boundaries.begin();
  const DD::Device::const_boundary_iterator end = boundaries.end();
  const map<const ElectricalContact*,
        double>::const_iterator bc_end = _simulation_voltages.end();
  for ( ; it != end; ++it)
  {
    if (_simulation_voltages.find(*it) == bc_end)
      _simulation_voltages[*it] = 0.0;

    _old_sim_voltages[*it] = 0.0;
  }

  _rebuild_eq_system = false;
}


void
DriftDiffusion::solve_newton(void) throw (PetscRuntimeError)
{

  PerfLog perf_log("solve_newton", false);

  // aliases for nicer code
  Options& params = get_options();
  SolverParameters& solver_params = params.solver_params;
  
  DD::Device& device = *_device;
  Mesh& mesh = _device->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  EquationSystems& equation_systems = get_equation_system();

  NonlinearImplicitSystem& system =
    equation_systems.get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

  NumericVector<Number>& solution = *(system.solution);
  NumericVector<Number>& old_solution =
    system.get_vector("old solution");


  // in 1D kspbcgs seems to work better than kspcgsl
  if (dim == 1)
    if (solver_params.ksp_type == KSPBCGSL)
      solver_params.ksp_type = KSPBCGS;


  // set the solver parameters (they could have changed since we made
  // the last calculation)
  set_solver_params(*system.nonlinear_solver);


  // set the right assembly function
  switch (_options.coupling)
  {
    case (POISSON | ECURRENT):
      system.nonlinear_solver->matvec = assemble<POISSON | ECURRENT>;
      break;
    case (POISSON | HCURRENT):
      system.nonlinear_solver->matvec = assemble<POISSON | HCURRENT>;
      break;
    case (CURRENTS):
      system.nonlinear_solver->matvec = assemble<CURRENTS>;
      break;
    case (POISSON):
      system.nonlinear_solver->matvec = assemble<POISSON>;
      break;
    case (ECURRENT):
      system.nonlinear_solver->matvec = assemble<ECURRENT>;
      break;
    case (HCURRENT):
      system.nonlinear_solver->matvec = assemble<HCURRENT>;
      break;
    default:
      system.nonlinear_solver->matvec = assemble<FULLYCOUPLED>;
  }


  //
  // solve for the desired contact voltages
  //
  ContactData voltages = _simulation_voltages;

  _options.C0_e = _options.n_max;
  _options.C0_h = _options.p_max;
  build_scaling();

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

      // if it did not converge, try half of the step
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



  // do mesh refinement if desired
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

      build_scaling();
      try { system.solve(); }
      catch (PetscDivergedError& e)
      {
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
DriftDiffusion::do_gummel_iterations(int max_it)
  throw (PetscRuntimeError, KSPDivergedError, SNESDivergedError)
{
  NonlinearImplicitSystem& system =
    get_equation_system().get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");
  
  SolverParameters& solver_params = get_options().solver_params;
  int nonlin_max_it = solver_params.nonlinear_max_iterations;

  try
  {
    for (int i = 0; i < max_it; i++)
    {
      system.nonlinear_solver->matvec = assemble<POISSON>;
      system.solve();

      _options.linearize_continuity_eq = true;
      solver_params.nonlinear_max_iterations = 1;
      set_solver_params(*system.nonlinear_solver);

      system.nonlinear_solver->matvec = assemble<ECURRENT>;
      system.solve();

      system.nonlinear_solver->matvec = assemble<HCURRENT>;
      system.solve();

      _options.linearize_continuity_eq = false;
      solver_params.nonlinear_max_iterations = nonlin_max_it;
      set_solver_params(*system.nonlinear_solver);
    }
  }
  catch (PetscRuntimeError err)
  {
    _options.linearize_continuity_eq = false;
    solver_params.nonlinear_max_iterations = nonlin_max_it;
    set_solver_params(*system.nonlinear_solver);
    throw(err);
  }
}

void
DriftDiffusion::solve_gummel(void) throw (PetscRuntimeError)
{
  do_gummel_iterations(get_options().max_gummel_iterations);
}


template <typename T>
void
DriftDiffusion::get_solution(const Elem* elem, const vector<Point>& p,
    vector<T>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  // this will contain the element in which p lies and for which
  // DriftDiffusion knows the potential
  const Elem* el = elem;
  
  // check if elem is in _element_list
  set<const Elem*>::iterator end = _element_list.end();
  set<const Elem*>::iterator it = _element_list.find(elem);

  if (it == end)
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();
    while (parent != NULL)
    {
      it = _element_list.find(parent);

      if (it != end)
        break; // we have found it, so get out of the while loop

      parent = parent->parent();
    }
    el = parent; // is NULL if no parent
  }

  if (el != NULL) // we found it!
    get_solution_secure(el, p, solution);
  else
  {
    // no parent, so check for children
    vector<const Elem*> tree;
    elem->family_tree(tree, false);

    set<const Elem*> elem_list;
    unsigned int len = tree.size();
    for (unsigned int i = 0; i < len; i++)
    {
      it = _element_list.find(tree[i]);
      if (it != end)
        elem_list.insert(tree[i]);
    }
    for (unsigned int i = 0; i < np; i++)
    {
      set<const Elem*>::iterator el_it = elem_list.begin();
      set<const Elem*>::iterator el_end = elem_list.end();
      for ( ; el_it != el_end; ++el_it)
      {
        el = *el_it;
        if (el->contains_point(p[i]))
        {
          get_solution_secure(el, p[i], solution[i]);
          // we have found it, so get out of the for loop
          break;
        }
      }
    }
  }
}


template <typename T>
void
DriftDiffusion::get_solution(const Elem* elem, const Point& p,
    T& solution)
{

  // this will contain the element in which p lies and for which
  // DriftDiffusion knows the potential
  const Elem* el = elem;
  
  // check if elem is in _element_list
  set<const Elem*>::iterator end = _element_list.end();
  set<const Elem*>::iterator it = _element_list.find(elem);

  if (it == end)
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();
    while (parent != NULL)
    {
      it = _element_list.find(parent);

      if (it != end)
        break; // we have found it, so get out of the while loop

      parent = parent->parent();
    }
    el = parent; // is NULL if no parent

    // no parent, so check for children
    if (el == NULL)
    {
      vector<const Elem*> tree;
      elem->family_tree(tree, false);
      
      unsigned int len = tree.size();
      for (unsigned int i = 0; i < len; i++)
      {
        it = _element_list.find(tree[i]);
        if (it != end)
        {
          if (tree[i]->contains_point(p))
          {
            // we have found it, so get out of the for loop
            el = tree[i];
            break;
          }
        }
      }
    }
  }
  // now el points to a valid element containing p or is NULL

  if (el != NULL)
    get_solution_secure(el, p, solution);
}

void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<double>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  NonlinearImplicitSystem* system;
  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  const NumericVector<Number>& ddsol = *(system->solution);

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");

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
      u  += phi[i][n] * ddsol(dof_indices_u[i]);

    // scale the potential back
    u  *= phi0;

    solution[n] = u;
  }
}

void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<DriftDiffusion::Solution>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  NonlinearImplicitSystem* system;
  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  const NumericVector<Number>& ddsol = *(system->solution);

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  for (unsigned int n = 0; n < np; n++)
  {
    double u = 0;
    double en = 0;
    double ep = 0;
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u  += phi[i][n] * ddsol(dof_indices_u[i]);
      en += phi[i][n] * ddsol(dof_indices_en[i]);
      ep += phi[i][n] * ddsol(dof_indices_ep[i]);
    }

    // scale the potential back
    u  *= phi0;
    en *= phi0;
    ep *= phi0;

    solution[n].potential = u;
    solution[n].fermi_e = en;
    solution[n].fermi_h = ep;
  }
}


void
DriftDiffusion::get_solution(const Elem* elem,
    std::vector<DriftDiffusion::Solution>& solution)
{
  NonlinearImplicitSystem* system =
    &_eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

  const DofMap& dof_map = system->get_dof_map();
  const NumericVector<Number>& sol = *(system->solution);

  const double phi0 = get_scaling().get_potential_scaling();
  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  assert(elem != NULL);

  solution.resize(elem->n_nodes());

  // get DOF indices
  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  for (unsigned int n = 0; n < elem->n_nodes(); n++)
  {
    solution[n].potential = phi0 * sol(dof_indices_u[n]);
    solution[n].fermi_e = phi0 * sol(dof_indices_en[n]);
    solution[n].fermi_h = phi0 * sol(dof_indices_ep[n]);
  }

}

/*
void
DriftDiffusion::get_electric_potential(const Elem* elem, const vector<Point>& p,
    vector<double>& potential)
{
  unsigned int np = p.size();
  potential.resize(np);
  if (np == 0) return;

}

double
DriftDiffusion::get_electric_potential(const Elem* elem, const Point& p)
{
  // this will contain the element in which p lies and for which
  // DriftDiffusion knows the potential
  const Elem* el = elem;
  
  // check if elem is in _element_list
  set<const Elem*>::iterator end = _element_list.end();
  set<const Elem*>::iterator it = _element_list.find(elem);

  if (it == end)
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();
    while (parent != NULL)
    {
      it = _element_list.find(parent);

      if (it != end)
        break; // we have found it, so get out of the while loop

      parent = parent->parent();
    }
    el = parent; // is NULL if no parent

    // no parent, so check for children
    if (el == NULL)
    {
      vector<const Elem*> tree;
      elem->family_tree(tree, false);
      
      unsigned int len = tree.size();
      for (unsigned int i = 0; i < len; i++)
      {
        it = _element_list.find(tree[i]);
        if (it != end)
        {
          if (tree[i]->contains_point(p))
          {
            // we have found it, so get out of the for loop
            el = tree[i];
            break;
          }
        }
      }
    }
  }
  // now el points to a valid element containing p or is NULL

  double u = 0;

  if (el != NULL)
  {
    // we found an element

    NonlinearImplicitSystem* system;
    system = &_eq_system->get_system<NonlinearImplicitSystem>(
        "drift-diffusion coupled");

    const NumericVector<Number>& solution = *(system->solution);

    const unsigned int dim = get_mesh().mesh_dimension();

    const DofMap& dof_map = system->get_dof_map();

    const unsigned int u_var = system->variable_number("potential");

    FEType fe_type = system->variable_type(u_var);
    AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

    vector<Point> point(1, FEInterface::inverse_map(dim, fe_type, el, p, 1e-6));
    fe->reinit(el, &point);

    vector<unsigned int> dof_indices_u;

    // element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    const unsigned int n_dofs = dof_indices_u.size();

    // the scaling parameters to scale back the result
    double phi0 = get_scaling().get_potential_scaling();

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
      u += phi[i][0] * solution(dof_indices_u[i]);

    // scale the potential back
    u *= phi0;
  }

  return u;
}*/


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
  
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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

    DriftDiffusionProperties* sc =
      device.get_element_data().get_data((*el)->top_parent());
    assert(sc != NULL);

    sc->reinit(elem);

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
          sc->calculate_all(u * phi0, en * phi0, ep * phi0, q_point[qp]);
        
          Real cond_e = Constants::e * sc->get_electron_conductivity();
          Real cond_h = Constants::e * sc->get_hole_conductivity();

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

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
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
  
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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

    DriftDiffusionProperties* sc =
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
        ElectricalContact* boundary_desc =
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
            sc->calculate_all(u * phi0, en * phi0, ep * phi0, q_point[qp]);
        
            Real cond_e = Constants::e * sc->get_electron_conductivity();
            Real cond_h = Constants::e * sc->get_hole_conductivity();

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
          
          // what is the outer normal in this point??
          // Idea: if x(s) > x(centroid), normal is +1
          //       else it is -1
          double x_c = elem->centroid()(0);
          double x_s = elem->point(s)(0);
          if (x_s < x_c)
          {
            dEfn = -dEfn;
            dEfp = -dEfp;
          }

          // calculate densities etc.
          sc->calculate_all(u * phi0, en * phi0, ep * phi0, elem->point(s));

          Real cond_e = Constants::e * sc->get_electron_conductivity();
          Real cond_h = Constants::e * sc->get_hole_conductivity();

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

void
DriftDiffusion::build_scaling(void)
{
  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);
  NumericVector<Number>& scaling = system->get_vector("scaling");

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_nodes();

  vector<unsigned short int> node_conn(nn);

  vector<double> nodal_val;

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  scaling.zero();

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

  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

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

      DriftDiffusionProperties* sc =
        device.get_element_data().get_data(elem->top_parent());
      assert(sc != NULL); 

      sc->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        sc->calculate_all(u, en, ep, elem->point(n));

        assert (node_conn[elem->node(n)] != 0);

        double nn  = sc->get_electron_density();
        double pp  = sc->get_hole_density();
        double ni = sc->get_intrinsic_density();
        double nn0 = sc->get_equilibrium_electron_density();
        double pp0 = sc->get_equilibrium_hole_density();
        double mue = sc->get_electron_mobility();
        double muh = sc->get_hole_mobility();
        double a = nn0 * exp(-nn / nn0 * 1);
        double b = pp0 * exp(-pp / pp0 * 1);
      
        scaling.set(dof_indices_u[n], 1);

        //double nodal_val = nn + a;
        double nodal_val = nn;
        nodal_val /= static_cast<Real>(node_conn[elem->node(n)]);
        scaling.add(dof_indices_en[n], nodal_val);

        //nodal_val = pp + b;
        nodal_val = pp;
        nodal_val /= static_cast<Real>(node_conn[elem->node(n)]);
        scaling.add(dof_indices_ep[n], nodal_val);
      }
    }
  }
}

void
DriftDiffusion::build_solution_vector(vector<double>& solution_vector)
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
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();

  const unsigned int n_vars  = 3;

  solution_vector.resize(nn * n_vars);
  vector<double> local(solution_vector.size());

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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

      assert(elem->n_nodes() == dof_indices_u.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        
        unsigned int id = n_vars * elem->node(n);
        local[id] = u;
        local[id + 1] = en;
        local[id + 2] = ep;
      }

    }
  }

#ifdef HAVE_MPI
  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  MPI_Allreduce (&local[0], &solution_vector[0], solution_vector.size(),
		 MPI_REAL, MPI_SUM, libMesh::COMM_WORLD);
#else
  solution_vector = local;
#endif
}

/*
void
DriftDiffusion::build_recombinations(vector<double>& densities,
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
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();

  const unsigned int n_vars  = 5;
  names.resize(n_vars);
  names[0] = "e_density";
  names[1] = "h_density";
  names[2] = "SRH";
  names[3] = "direct";
  names[4] = "Auger";

  vector<double> recomb(3, 0.0);

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

      DriftDiffusionProperties* sc =
        device.get_element_data().get_data(elem->top_parent());
      assert(sc != NULL); 

      sc->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      //sol.resize(dof_indices.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        
        sc->calculate_all(u, en, ep, elem->point(n));
        sc->get_net_recombination_rates(recomb);

        assert (node_conn[elem->node(n)] != 0);

        unsigned int id = n_vars * elem->node(n);
        double nodal_val = sc->get_electron_density();
        local[id] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = sc->get_hole_density();
        local[id + 1] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        //nodal_val = sc->get_net_electron_recombination_rate();
        nodal_val = recomb[0];
        local[id + 2] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);
        
        nodal_val = recomb[1];
        local[id + 3] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = recomb[2];
        local[id + 4] +=
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
*/


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
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();

  const unsigned int n_vars  = 8;
  names.resize(n_vars);
  names[0] = "electron_density";
  names[1] = "hole_density";
  names[2] = "ionized_donors";
  names[3] = "ionized_acceptors";
  names[4] = "total_charge";
  names[5] = "SRH";
  names[6] = "direct";
  names[7] = "Auger";

  vector<double> recomb(3, 0.0);

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

      DriftDiffusionProperties* sc =
        device.get_element_data().get_data(elem->top_parent());
      assert(sc != NULL); 

      sc->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      //sol.resize(dof_indices.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        
        sc->calculate_all(u, en, ep, elem->point(n));
        sc->get_net_recombination_rates(recomb);

        assert (node_conn[elem->node(n)] != 0);

        unsigned int id = n_vars * elem->node(n);
        double nodal_val = sc->get_electron_density();
        local[id] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = sc->get_hole_density();
        local[id + 1] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = sc->get_ionized_donor_density();
        local[id + 2] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);
        
        nodal_val = sc->get_ionized_acceptor_density();
        local[id + 3] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = sc->get_charge_density();
        local[id + 4] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = recomb[0];
        local[id + 5] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);
        
        nodal_val = recomb[1];
        local[id + 6] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = recomb[2];
        local[id + 7] +=
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
DriftDiffusion::build_current_density(vector<double>& current,
    vector<string>& names)
{
  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;
  
  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_elem();

  const unsigned int n_vars  = 10;
  names.resize(n_vars);
  names[0] = "Jn_x";
  names[1] = "Jn_y";
  names[2] = "Jn_z";
  names[3] = "Jp_x";
  names[4] = "Jp_y";
  names[5] = "Jp_z";
  names[6] = "J_x";
  names[7] = "J_y";
  names[8] = "J_z";
  names[9] = "J";

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
    current[id]     = phi0 / x0 * sigma_e * en_x;
    current[id + 1] = phi0 / x0 * sigma_e * en_y;
    current[id + 2] = phi0 / x0 * sigma_e * en_z;
    current[id + 3] = phi0 / x0 * sigma_h * ep_x;
    current[id + 4] = phi0 / x0 * sigma_h * ep_y;
    current[id + 5] = phi0 / x0 * sigma_h * ep_z;
    double jx = phi0 / x0 * (sigma_e * en_x + sigma_h * ep_x);
    double jy = phi0 / x0 * (sigma_e * en_y + sigma_h * ep_y);
    double jz = phi0 / x0 * (sigma_e * en_z + sigma_h * ep_z);
    current[id + 6] = jx;
    current[id + 7] = jy;
    current[id + 8] = jz;
    current[id + 9] = sqrt(jx * jx + jy * jy + jz * jz);

    elem_number++;
  }
  current.resize(elem_number * n_vars);
}


void
DriftDiffusion::build_electric_field(vector<double>& field,
    vector<string>& names)
{
  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_elem();

  const unsigned int n_vars  = 3;
  names.resize(n_vars);
  names[0] = "electric_field_x";
  names[1] = "electric_field_y";
  names[2] = "electric_field_z";

  field.resize(nn * n_vars);

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();
  const double x0 = get_options().mesh_units;

  const unsigned int u_var = system->variable_number("potential");
  
  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, libMeshEnums::CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  vector<unsigned int> dof_indices_u;

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    DriftDiffusionProperties* sc =
      device.get_element_data().get_data(elem->top_parent());
    assert(sc != NULL); 

    sc->reinit(elem);

    //vector<Point> centroid(1, elem->centroid());
    fe->reinit(elem);
    
    unsigned int n_dofs = dof_indices_u.size();
    // get the solution values at the centroid
    Real Ex = 0.0;
    Real Ey = 0.0;
    Real Ez = 0.0;
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      Ex  += dphi[i][0](0) * solution(dof_indices_u[i]);
      Ey  += dphi[i][0](1) * solution(dof_indices_u[i]);
      Ez  += dphi[i][0](2) * solution(dof_indices_u[i]);
    }

    unsigned int id = n_vars * elem_number;
    field[id] = -phi0 * Ex / x0;
    field[id + 1] = -phi0 * Ey / x0;
    field[id + 2] = -phi0 * Ez / x0;

    elem_number++;
  }
  field.resize(elem_number * n_vars);
}


void
DriftDiffusion::build_elem_band_edges(vector<double>& field,
    vector<string>& names)
{
  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
  NonlinearImplicitSystem* system;

  system = &_eq_system->get_system<NonlinearImplicitSystem>(
      "drift-diffusion coupled");

  // aliases for nicer code
  const DD::Device& device = *(_device);
  const Mesh& mesh = _device->get_mesh();
  const NumericVector<Number>& solution = *(system->solution);

  //const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  unsigned int nn  = 0;

  MeshBase::const_element_iterator it1 =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end1 =
    mesh.active_elements_end(); 
  for ( ; it1 != end1; ++it1)   {nn++;}

  const unsigned int n_vars  = 3;
  names.resize(n_vars);
  names[0] = "conduction_band_edge";
  names[1] = "valence_band_edge";
  names[2] = "equilibrium_fermi_level";

  field.resize(nn * n_vars);

  // the scaling parameters to scale back the result
  // double phi0 = get_scaling().get_potential_scaling();
  //const double x0 = get_options().mesh_units;

  // const unsigned int u_var = system->variable_number("potential");
  
  //FEType fe_type = system->variable_type(u_var);
  //AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  // QGauss qrule(dim, libMeshEnums::CONSTANT);
  //fe->attach_quadrature_rule(&qrule);

  //vector<unsigned int> dof_indices_u;

  // element shape function gradients
  //const vector<vector<Real> >& phi = fe->get_phi();

  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    //dof_map.dof_indices(elem, dof_indices_u, u_var);

    DriftDiffusionProperties* sc =
      device.get_element_data().get_data(elem->top_parent());

    assert(sc != NULL); 

    sc->reinit(elem);

    //    fe->reinit(elem);

    
    
    // unsigned int n_dofs = dof_indices_u.size();
    // get the solution values at the centroid
    // Real u = 0.0;
    // for (unsigned int i = 0; i < n_dofs; i++)
    //  u  += phi[i][0] * solution(dof_indices_u[i]);

    unsigned int id = n_vars * elem_number;
    //field[id] = sc->get_conduction_band_edge() - phi0 * u;
    //field[id + 1] = sc->get_valence_band_edge() - phi0 * u;
    field[id] = sc->get_conduction_band_edge();
    field[id + 1] = sc->get_valence_band_edge();
    field[id + 2] = sc->get_equilibrium_fermi_level();

    elem_number++;
  }
  //field.resize(elem_number * n_vars);
}


void
DriftDiffusion::build_band_edges(vector<double>& band_edges,
    vector<string>& names)
{
  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
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

  const unsigned int n_vars  = 5;
  names.resize(n_vars);
  names[0] = "conduction_band_edge";
  names[1] = "valence_band_edge";
  names[2] = "electron_quasi_fermi_level";
  names[3] = "hole_quasi_fermi_level";
  names[4] = "electrical_potential";

  band_edges.resize(nn * n_vars);

  vector<double> local(band_edges.size());
  vector<unsigned short int> node_conn(nn);

  vector<double> nodal_val;

  fill(band_edges.begin(), band_edges.end(), 0.0);
  fill(local.begin(), local.end(), 0.0);

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

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

      DriftDiffusionProperties* sc =
        device.get_element_data().get_data(elem->top_parent());
      assert(sc != NULL); 

      sc->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      //sol.resize(dof_indices.size());

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Real u  = phi0 * solution(dof_indices_u[n]);
        Real en = phi0 * solution(dof_indices_en[n]);
        Real ep = phi0 * solution(dof_indices_ep[n]);
        
        sc->calculate_all(u, en, ep, elem->point(n));

        assert (node_conn[elem->node(n)] != 0);

        unsigned int id = n_vars * elem->node(n);
        double nodal_val = sc->get_conduction_band_edge() - u;
        local[id] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        nodal_val = sc->get_valence_band_edge() - u;
        local[id + 1] +=
          nodal_val / static_cast<Real>(node_conn[elem->node(n)]);

        local[id + 2] = en;
        local[id + 3] = ep;
        local[id + 4] = u;
      }
    }
  }

#ifdef HAVE_MPI
  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  MPI_Allreduce (&local[0], &band_edges[0], band_edges.size(),
		 MPI_REAL, MPI_SUM, libMesh::COMM_WORLD);
#else
  densities = local;
#endif
}


// TODO
// This routine needs some optimization, mostly concerning the
// scaling
template <int coupling>
void
DriftDiffusion::assemble(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  PerfLog perf_log("Matrix assembly", false);
  perf_log.start_event("assembly");
  
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
  bool linearize = options.linearize_continuity_eq;

  ContactData& simulation_voltages = _this->_simulation_voltages;
  BoundaryNodeList& dirichlet_nodes = _this->_dirichlet_nodes;


  //
  // some scaling stuff...
  // 
  // NOTE: the mesh and all paramters were not explicitly scaled, so
  //       we have to treat scaling by explicit division/multiplication
  //       
  // maximum density of electrons
  double n_max = 1;
  // maximum density of holes
  double p_max = 1;
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
  // density scaling for electrons
  double C0_e = options.C0_e;
  // density scaling for holes
  double C0_h = options.C0_h;
  // scaling for recombination rates
  double R0_e = C0_e / scaling.get_time_scaling();
  double R0_h = C0_h / scaling.get_time_scaling();
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

    sc->reinit(elem);

    // get the nodal scaling values
    vector<double> n0(n_dofs, 1);
    vector<double> p0(n_dofs, 1);
    if (options.local_scaling)
    {
      NumericVector<Number>& scaling = system.get_vector("scaling");
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        n0[i] = scaling(dof_indices_en[i]) / C0_e;
        p0[i] = scaling(dof_indices_ep[i]) / C0_h;
      }
    }

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
      sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, q_point[qp]);
      double n = sc->get_electron_density();
      double p = sc->get_hole_density();
      double epsilon = sc->get_relative_permittivity();
      double l2_eps = l2 * epsilon;

      double Rn = sc->get_net_electron_recombination_rate();
      //Rn = (fabs(Rn) < 1.0e-19) ? 0.0 : Rn;
      double Rp = sc->get_net_hole_recombination_rate();
      //Rp = (fabs(Rp) < 1.0e-19) ? 0.0 : Rp;
      
      // remember the maximum densities
      n_max = (n_max > n) ? n_max : n;
      p_max = (p_max > p) ? p_max : p;

      double ni = sc->get_intrinsic_density();
      double nn0 = sc->get_equilibrium_electron_density();
      double pp0 = sc->get_equilibrium_hole_density();
      double mue = sc->get_electron_mobility();
      double muh = sc->get_hole_mobility();
      double a = 0.0;
      double b = 0.0;
      if (options.artificial_drift)
      {
        a = exp(-n / nn0 * 1);
        b = exp(-p / pp0 * 1);
      }

      // NOTE: sigma_e = mu_e * n is the electron conductivity
      double sigma_e = (sc->get_electron_conductivity()
          + mue * nn0 * a) / (mu0 * C0_e);
      double sigma_h = (sc->get_hole_conductivity()
          + muh * pp0 * b) / (mu0 * C0_h);
    
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
            J * (dphi[i][qp] * dphi[j][qp]) * x0_mesh * x0_mesh;
          
          if (coupling & POISSON)
            Kuu(i,j) += l2_eps * laplace;
          
          if (coupling & ECURRENT)
            Knn(i,j) += sigma_e * laplace / n0[i];
          
          if (coupling & HCURRENT)
            Kpp(i,j) += sigma_h * laplace / p0[i];
        }

        if (!(coupling & POISSON))
          Kuu(i,i) += 1;
        
        if (!(coupling & ECURRENT))
          Knn(i,i) += 1;
        
        if (!(coupling & HCURRENT))
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
          //dRn[id] = phi0 / R0_e
          dRn[id] = phi0 / R0_e * (1 - a)
            * sc->get_net_electron_recombination_rate_derivatives()[id];
          //dRp[id] = phi0 / R0_h
          dRp[id] = phi0 / R0_h * (1 - b)
            * sc->get_net_hole_recombination_rate_derivatives()[id];
        }
        //if (fabs(Rn) < 1e-3)
        if (Rn == 0.0)
          dRn[0] = dRn[1] = dRn[2] = 0.0;
        //if (fabs(Rp) < 1e-3)
        if (Rp == 0.0)
          dRp[0] = dRp[1] = dRp[2] = 0.0;

        // d(sigma_n)/du * element-jacobian
        // sigma_n = mu_n * n means the conductivity of electrons
        Real dsigma_e = J * phi0 / (mu0 * C0_e)
          * sc->get_electron_conductivity_derivatives()[0]
               * (1 - a);
        Real dsigma_h = J * phi0 / (mu0 * C0_h)
          * sc->get_hole_conductivity_derivatives()[0]
               * (1 - b);

        if (linearize)
        {
          dRn[0] = dRn[1] = dRn[2] = 0.0;
          dRp[0] = dRp[1] = dRp[2] = 0.0;
          dsigma_e = dsigma_h = 0.0;
        }

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            
            Real dsigma_e_x_phi = dsigma_e * phi[j][qp] / n0[i];
            Real dsigma_h_x_phi = dsigma_h * phi[j][qp] / p0[i];
            for (unsigned int k = 0; k < n_dofs; k++)
            {
              Real laplace = (dphi[i][qp] * dphi[k][qp]) * x0_mesh * x0_mesh;
            
              if (coupling & ECURRENT)
              {
                Real elem_contrib =
                  dsigma_e_x_phi * laplace * Xn(k);

                if (coupling & POISSON)
                  Knu(i,j) += elem_contrib;

                Knn(i,j) += elem_contrib;
              }

              if (coupling & HCURRENT)
              {
                Real elem_contrib =
                  dsigma_h_x_phi * laplace * Xp(k);
                
                if (coupling & POISSON)
                  Kpu(i,j) += elem_contrib;

                Kpp(i,j) += elem_contrib;
              }
            }

            // The dFe_i/dX_j part
            Real phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (coupling & POISSON)
            {
              Kuu(i,j) -= drho[0] * phi_i_x_phi_j;
              
              if (coupling & ECURRENT)
                Kun(i,j) -= drho[1] * phi_i_x_phi_j;

              if (coupling & HCURRENT)
                Kup(i,j) -= drho[2] * phi_i_x_phi_j;
            }            
            
            if (coupling & ECURRENT)
            {
              if (coupling & POISSON)
                Knu(i,j) += dRn[0] * phi_i_x_phi_j / n0[i];

              Knn(i,j) += dRn[1] * phi_i_x_phi_j / n0[i];

              if (coupling & HCURRENT)
                Knp(i,j) += dRn[2] * phi_i_x_phi_j / n0[i];
            }

            if (coupling & HCURRENT)
            {
              if (coupling & POISSON)
                Kpu(i,j) -= dRp[0] * phi_i_x_phi_j / p0[i];

              if (coupling & ECURRENT)
                Kpn(i,j) -= dRp[1] * phi_i_x_phi_j / p0[i];
              
              Kpp(i,j) -= dRp[2] * phi_i_x_phi_j / p0[i];
            }

          }
        }
      }

      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        Real J_x_rho = J * sc->get_charge_density() / C0;
        Real J_x_P0 = J / P0;

        // net recombination rate
        Real J_x_Rn = J * Rn / R0_e * (1 - a);
        Real J_x_Rp = J * Rp / R0_h * (1 - b);

        RealVectorValue P(sc->get_total_polarization());
        P *= J_x_P0;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          Real net_recomb_e = J_x_Rn * phi[i][qp] / n0[i];
          Real net_recomb_h = J_x_Rp * phi[i][qp] / p0[i];
          
          if (coupling & POISSON)
            Fu(i) -= J_x_rho * phi[i][qp] + (P * dphi[i][qp]) * x0_mesh;
          else
            Fu(i) -= Xu(i);
          
          if (coupling & ECURRENT)
            Fn(i) += net_recomb_e;
          else
            Fn(i) -= Xn(i);

          if (coupling & HCURRENT)
            Fp(i) -= net_recomb_h;
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
        ElectricalContact* contact = 
          device.get_boundary_data().get_data(side);

        // for von Neumann or mixed type boundary conditions
        vector<double> coeff(3, 0.0);
        vector<double> value(3, 0.0);

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
            sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, q_point[qp]);
            double epsilon = sc->get_relative_permittivity();
            double l2_eps = l2 * epsilon;

            // get the boundary condition coefficients
            if (contact != NULL)
            {
              contact->set_material(sc);
              double a, c;

              if (coupling & POISSON)
              {
                contact->get_normal_derivative(POTENTIAL, a, c);
                coeff[0] = a * x0;
                value[0] = c * x0 / phi0;
              }
              if (coupling & ECURRENT)
              {
                contact->get_normal_derivative(FERMIE, a, c);
                coeff[1] = a * x0;
                value[1] = c * x0 / phi0;
              }
              if (coupling & HCURRENT)
              {
                contact->get_normal_derivative(FERMIH, a, c);
                coeff[2] = a * x0;
                value[2] = c * x0 / phi0;
              }
            }



            // the jacobian x weight x scaling
            double J = JxW_face[qp] / Jface_scale;

            // first the contributions to Ke_ij
            for (unsigned int i = 0; i < n_dofs; i++)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {

                Real phi_i_x_phi_j =
                  J * phi_face[i][qp] * phi_face[j][qp];

                if (coupling & POISSON)
                  Kuu(i,j) += l2_eps * coeff[0] * phi_i_x_phi_j;

                if (coupling & ECURRENT)
                  Knn(i,j) += coeff[1] * phi_i_x_phi_j;

                if (coupling & HCURRENT)
                  Kpp(i,j) += coeff[2] * phi_i_x_phi_j;
              }
            }

            // contribution to the jacobian
            if (jacobian != NULL)
            {
            }

            // contribution to -Fe_i
            if (residual != NULL)
            {
              RealVectorValue P(sc->get_total_polarization());
              double Pn = (P * face_normals[qp]) / P0;
              double value_u = J * (l2_eps * value[0] - Pn);
              double value_n = J * value[1] / (mu0 * C0_e);
              double value_p = J * value[2] / (mu0 * C0_h);

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                if (coupling & POISSON)
                  Fu(i) -= value_u * phi_face[i][qp];

                if (coupling & ECURRENT)
                  Fn(i) -= value_n * phi_face[i][qp] / n0[i];

                if (coupling & HCURRENT)
                  Fp(i) -= value_p * phi_face[i][qp] / p0[i];
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
          sc->calculate_all(phi0 * u, phi0 * en, phi0 * ep, elem->point(s));
          double epsilon = sc->get_relative_permittivity();
          double l2_eps = l2 * epsilon;

          // get the boundary condition coefficients
          if (contact != NULL)
          {
            contact->set_material(sc);
            double a, c;

            if (coupling & POISSON)
            {
              contact->get_normal_derivative(POTENTIAL, a, c);
              coeff[0] = a * x0;
              value[0] = c * x0 / phi0;
            }
            if (coupling & ECURRENT)
            {
              contact->get_normal_derivative(FERMIE, a, c);
              coeff[1] = a * x0;
              value[1] = c * x0 / phi0;
            }
            if (coupling & HCURRENT)
            {
              contact->get_normal_derivative(FERMIH, a, c);
              coeff[2] = a * x0;
              value[2] = c * x0 / phi0;
            }
          }


          // first the contributions to Ke_ij
          if (coupling & POISSON)
            Kuu(s,s) += l2_eps * coeff[0];

          if (coupling & ECURRENT)
            Knn(s,s) += coeff[1];

          if (coupling & HCURRENT)
            Kpp(s,s) += coeff[2];

          // contribution to -Fe_i
          if (residual != NULL)
          {
            double Pn =  sc->get_total_polarization()(0) / P0;
            // what is the outer normal in this point??
            // Idea: if x(s) > x(centroid), normal is +1
            //       else it is -1
            double x_c = elem->centroid()(0);
            double x_s = elem->point(s)(0);
            Pn = (x_s > x_c) ? Pn : -Pn;
            double value_u = l2_eps * value[0] - Pn;
            double value_n = value[1] / (mu0 * C0_e * n0[s]);
            double value_p = value[2] / (mu0 * C0_h * p0[s]);

            if (coupling & POISSON)
              Fu(s) -= value_u;

            if (coupling & ECURRENT)
              Fn(s) -= value_n;

            if (coupling & HCURRENT)
              Fp(s) -= value_p;
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
          ElectricalContact* contact = node_it->second;
          contact->set_material(sc);

          if (coupling & POISSON)
          {
            if (contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
            {
              double val = (contact->get_boundary_value(POTENTIAL)
                  + simulation_voltages[contact]) / phi0;
              Ke.condense(i, i, -val, Fe);
            }
            else if (contact->get_type(POTENTIAL) == ElectricalContact::PINNING)
            {
              double val = contact->get_boundary_value(POTENTIAL) / phi0;
              Ke.condense(i, i, -val, Fe);
              Ke(i, i + n_dofs) = 1.0;
            }
          }

          if (coupling & ECURRENT)
          {
            if (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
            {
              double val = (contact->get_boundary_value(FERMIE)
                  - simulation_voltages[contact]) / phi0;
              Ke.condense(i + n_dofs, i + n_dofs, -val, Fe);
            }
          }

          if (coupling & HCURRENT)
          {
            if (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET)
            {
              double val = (contact->get_boundary_value(FERMIH)
                  - simulation_voltages[contact]) / phi0;
              Ke.condense(i + 2 * n_dofs, i + 2 * n_dofs, -val, Fe);
            }
          }
        }
      }
    }
    else
    {
      // TODO this needs to be checked!!!

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
          if (dof_map.is_constrained_dof(dof_indices_u[i]))
            is_done = false;

          node_it = dirichlet_nodes.find(parent->get_node(i));
          if (node_it != end)
          {
            ElectricalContact* contact = node_it->second;
            contact->set_material(sc);

            // loop over all DOFs occurring in the constrained matrix
            for (unsigned int id = 0; id < n_dofs_tot; id++)
            {

              if (coupling & POISSON)
              {
                if (contact->get_type(POTENTIAL) ==
                    ElectricalContact::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_u[i])
                  {
                    double val = (contact->get_boundary_value(POTENTIAL)
                      + simulation_voltages[contact]) / phi0;
                    Ke.condense(id, id, -val, Fe);
                  }
                }
                else if (contact->get_type(POTENTIAL) ==
                    ElectricalContact::PINNING)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_u[i])
                  {
                    double val = contact->get_boundary_value(POTENTIAL) / phi0;
                    Ke.condense(i, i, -val, Fe);
                    Ke(id, id + n_dofs) = 1.0;
                  }
                }
              }

              if (coupling & ECURRENT)
              {
                if (contact->get_type(FERMIE) ==
                    ElectricalContact::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_en[i])
                  {
                    double val = (contact->get_boundary_value(FERMIE)
                      - simulation_voltages[contact]) / phi0;
                    Ke.condense(id, id, -val, Fe);
                  }
                }
              }

              if (coupling & HCURRENT)
              {
                if (contact->get_type(FERMIH) ==
                    ElectricalContact::DIRICHLET)
                {
                  // is it a boundary DOF?
                  if (dof_indices[id] == dof_indices_ep[i])
                  {
                    double val = (contact->get_boundary_value(FERMIE)
                      - simulation_voltages[contact]) / phi0;
                    Ke.condense(id, id, -val, Fe);
                  }
                }
              }

            } // end loop over all DOFs 
          }
        } // end loop over the nodes of the parent element
      }
    }

    perf_log.start_event("add");
    if (residual != NULL)
    {
      for (unsigned int i = 0; i < n_dofs_tot; i++)
        for (unsigned int j = 0; j < n_dofs_tot; j++)
          Fe(i) += Ke(i,j) * x(dof_indices[j]);

      residual->add_vector(Fe, dof_indices);
    }
    else
      jacobian->add_matrix(Ke, dof_indices);
    perf_log.stop_event("add");

  } // end loop over elements

  // put the maximum densities back into the options
  options.n_max = n_max;
  options.p_max = p_max;
  
  perf_log.stop_event("assembly");
} 


//
// explicit instantiations of template methods
//
template void
DriftDiffusion::get_solution<double>(const Elem* elem, const Point& p,
    double& solution);

template void
DriftDiffusion::get_solution<double>(const Elem* elem, const vector<Point>& p,
    vector<double>& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::Solution>(const Elem* elem,
    const vector<Point>& p, vector<DriftDiffusion::Solution>& solution);
