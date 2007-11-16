// $Id$

// module includes
#include "DriftDiffusion.h"
#include "SimulationEnvironment.h"
#include "Scaling.h"
#include "Material.h"
#include "Boundary.h"
#include "ElectricalContact.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSolver.h"
#include "SolveFailedException.h"
#include "GraceIO.h"

// libmesh includes
#include "node.h"
#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "equation_systems.h"
#include "mesh_refinement.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

// C++ includes

using namespace std;
using namespace DriftDiffusionDefs;


//
// Module interface
//

TIBER_MODULE(DriftDiffusion,driftdiffusion)



DriftDiffusion*
DriftDiffusion::_this;



DriftDiffusion::Options::Options(void)
  : mesh_refinement(false),
    max_refinement_steps(5),
    max_refinement_level(8),
    refine_fraction(0.7),
    coarsen_fraction(0.3),
    refinement_tolerance(1e-6),
    integration_order(libMeshEnums::FIFTH),
    solver_method(NEWTON),
    max_gummel_iterations(5),
    scaling_type(Scaling::UNITS),
    coupling(FULLYCOUPLED),
    scheme(FEM),
    current_calculation(RSTF),
    exact_newton(true)
{
}





DriftDiffusion::Options::Options(const Options& rhs)
  : mesh_refinement(rhs.mesh_refinement),
    max_refinement_steps(rhs.max_refinement_steps),
    max_refinement_level(rhs.max_refinement_level),
    refine_fraction(rhs.refine_fraction),
    coarsen_fraction(rhs.coarsen_fraction),
    refinement_tolerance(rhs.refinement_tolerance),
    integration_order(rhs.integration_order),
    solver_method(rhs.solver_method),
    max_gummel_iterations(rhs.max_gummel_iterations),
    solver_params(rhs.solver_params),
    scaling_type(rhs.scaling_type),
    coupling(rhs.coupling),
    scheme(rhs.scheme),
    current_calculation(rhs.current_calculation),
    exact_newton(rhs.exact_newton)
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
    integration_order = rhs.integration_order;
    solver_method = rhs.solver_method;
    max_gummel_iterations = rhs.max_gummel_iterations;
    solver_params = rhs.solver_params;
    scaling_type = rhs.scaling_type;
    coupling = rhs.coupling;
    scheme = rhs.scheme;
    current_calculation = rhs.current_calculation;
    exact_newton = rhs.exact_newton;
  }
  return *this;
}





DriftDiffusion::SolverParameters::SolverParameters(void)
  : nonlinear_tolerance(1e-9), 
    nonlinear_abs_tolerance(1e-15),
    nonlinear_step_tolerance(1e-3),
    nonlinear_max_iterations(20),
    linear_tolerance(1e-6),
    linear_abs_tolerance(1e-50),
    linear_max_iterations(500),
    ls_maxstep(1.0),
    ls_type(3),
    ksp_type(BICGSTAB),
    pc_type(ILU_PRECOND),
    nonlinear_solver("petsc"),
    linear_solver("petsc")
{
  // TODO read default values from some text file
}





DriftDiffusion::SolverParameters::SolverParameters(const SolverParameters& rhs)
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
    pc_type(rhs.pc_type),
    nonlinear_solver(rhs.nonlinear_solver),
    linear_solver(rhs.linear_solver)
{
}





DriftDiffusion::SolverParameters&
DriftDiffusion::SolverParameters::operator=(const SolverParameters& rhs)
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
    nonlinear_solver = rhs.nonlinear_solver;
    linear_solver = rhs.linear_solver;
  }
  return *this;
}




DriftDiffusion::DriftDiffusion(void)
  : _rebuild_eq_system(true),
    _linear_solver(NULL)
{
}




DriftDiffusion::~DriftDiffusion(void)
{
  cleanup_solver();
}
 



PhysicalModel*
DriftDiffusion::create_physical_model(const ModelOptions& options) const
throw (ModelErrorException)
{
  const string& modelname = options.get_option("model", "unstrained");

  DriftDiffusionProperties* model =
    DriftDiffusionProperties::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "DriftDiffusion: No such physical model: " + modelname);

  return model;
}





BoundaryProperties*
DriftDiffusion::create_boundary_model(const ModelOptions& options) const
throw (ModelErrorException)
{
  const string& modelname = options.get_option("type", "");

  ElectricalContact* model = NULL;

  if (modelname != "")
  {
    model =  ElectricalContact::create(modelname, options);

    if (model == NULL)
      throw ModelErrorException(
          "DriftDiffusion: No such boundary model: " + modelname);
  }

  return model;
}




void
DriftDiffusion::compute_scaling(Scaling::ScalingType type)
{
  if (type == Scaling::NONE)
  {
    get_scaling().set_scaling_type(type);
    get_scaling().set_potential_scaling(1);
    get_scaling().set_length_scaling(1);
    get_scaling().set_mobility_scaling(1);
    get_scaling().set_density_scaling(1);
    return;
  }

  // we calculate in cm!
  double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  
  // the scaling parameters should never be zero
  // they are in any case positive, so it will
  // always find the maximum
  double x0 = -1;
  double phi0 = SimulationOptions::T * Constants::k_B;
  double mu0 = -1;
  double C0 = 1;
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

    assert(_device->get_material(elem->subdomain_id()) != NULL);
    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          _device->get_material(elem->subdomain_id())->get_model(get_id()));

    sc->reinit(elem);
    sc->set_coordinates(elem->centroid());
    sc->set_potentials(sc->get_equilibrium_fermi_level());
    sc->set_electric_field(RealGradient(0));

    sc->calculate_densities();
    sc->calculate_ionized_dopants();
    sc->calculate_mobilities();

    // TODO get max of polarisation
    
    double mu = sc->get_hole_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;
    mu = sc->get_electron_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;

    // I don't know what is better...
    double C = fabs(sc->get_material()->get_net_doping_density());
    //double C = fabs(sc->get_ionized_donor_density() -
    //    sc->get_ionized_acceptor_density());
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

  get_scaling().set_scaling_type(type);
  get_scaling().set_potential_scaling(phi0);
  //get_scaling.set_potential_scaling(1.0);
  get_scaling().set_length_scaling(x0 * mesh_units);
  get_scaling().set_mobility_scaling(mu0);
  get_scaling().set_density_scaling(C0);
}




void
DriftDiffusion::set_electron_fermi_level(double Ef_n)
{
  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  NumericVector<Number>& solution = system.get_solution_vector();
  
  const unsigned int var = system.variable_number("fermi_e");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = Ef_n / phi0;

  Mesh& mesh = get_mesh();
  Mesh::element_iterator it = mesh.active_elements_begin();
  const Mesh::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      unsigned int id = 
        elem->get_node(i)->dof_number(system.number(), var, 0);
      solution.set(id, level);
    }
  }
}




void
DriftDiffusion::set_hole_fermi_level(double Ef_p)
{
  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  NumericVector<Number>& solution = system.get_solution_vector();
  

  const unsigned int var = system.variable_number("fermi_h");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = Ef_p / phi0;

  Mesh& mesh = get_mesh();
  Mesh::element_iterator it = mesh.active_elements_begin();
  const Mesh::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      unsigned int id = 
        elem->get_node(i)->dof_number(system.number(), var, 0);
      solution.set(id, level);
    }
  }
}




void
DriftDiffusion::set_electric_potential(double pot)
{
  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  NumericVector<Number>& solution = system.get_solution_vector();
  
  const unsigned int var = system.variable_number("potential");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = -pot / phi0;

  Mesh& mesh = get_mesh();
  Mesh::element_iterator it = mesh.active_elements_begin();
  const Mesh::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      unsigned int id = 
        elem->get_node(i)->dof_number(system.number(), var, 0);
      solution.set(id, level);
    }
  }
}




void
DriftDiffusion::find_dirichlet_nodes(void)
{

  SimulationEnvironment& env = get_environment();

  SimulationEnvironment::BoundaryNodeIterator it(env.boundary_nodes_begin());
  SimulationEnvironment::BoundaryNodeIterator end(env.boundary_nodes_end());
  
  for ( ; it != end; ++it)
  {
    ID id = it->second;

    Boundary* bd = env.get_boundary(id);

    ElectricalContact* contact = NULL;
    if (bd != NULL)
      contact = dynamic_cast<ElectricalContact*>(
          bd->get_boundary_properties(get_id()));

    if (contact != NULL)
    {
      if ((contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
          || (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
          || (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET))
      {
        _dirichlet_nodes[it->first] = bd;
      }
    }
  }
}



void
DriftDiffusion::find_dielectric_boundary_nodes(void)
{
  Mesh& mesh = get_mesh();
  Mesh::element_iterator it = mesh.active_elements_begin();
  const Mesh::element_iterator end = mesh.active_elements_end();
  
  for ( ; it != end; ++it)
  {
    const Elem* el = *it;

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          _device->get_material(el->subdomain_id())->get_model(get_id()));
    
    // we are only interested in boundaries between semiconductor/dielectric
    if (sc->is_dielectric())
    {
      for (unsigned s = 0; s < el->n_sides(); s++)
      {
        if (get_environment().is_inner_boundary(ElementSide(el, s)))
        {
          // get the model of the neighbor element
          DriftDiffusionProperties* scn =
            dynamic_cast<DriftDiffusionProperties*>(
                _device->get_material(
                  el->neighbor(s)->subdomain_id())->get_model(get_id()));

          // if neighbor is not dielectric we record it
          if (!scn->is_dielectric())
          {
            AutoPtr<Elem> side(el->build_side(s));
            for (unsigned int i = 0; i < side->n_nodes(); i++)
              _dielectric_boundary_nodes.insert(side->get_node(i));
          }
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
    get_equation_systems().delete_system(get_equation_system_name());
    _rebuild_eq_system = true;
  }
}




void
DriftDiffusion::cleanup_solver(void)
{

  // erase boundary current data structure
  _boundary_currents.erase(_boundary_currents.begin(),
      _boundary_currents.end());

  // erase dirichlet nodes data structure
  _dirichlet_nodes.erase(_dirichlet_nodes.begin(),
      _dirichlet_nodes.end());

  delete _linear_solver;

  reset_solver();
}




void
DriftDiffusion::do_solve(void)
{
  
  cout << endl;
  cout << "<<-------------------------------------------------------------------"
    << endl;
  cout << "DriftDiffusion (name: " << get_name() << ")" << endl;

  // rebuild the system if needed
  //rebuild_equation_system();

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;
  
  parse_options();

  /* This is really stupid! How do you want to make selfconsistent loops ???
    // we do not solve anything if the simulation voltages are the
    // same as before
  */

  bool equilibrium = true;
  // this does not what we think in some cases!
  //bool accept_failure = true;
  bool same_potentials = true;
  bool accept_failure = false;
  {
    ContactData::iterator it(_voltages.begin());
    const ContactData::iterator end(_voltages.end());
    for ( ; it != end; ++it)
    {
      ElectricalContact* bd = static_cast<ElectricalContact*>(
          it->first->get_boundary_properties(get_id()));

      double voltage = bd->get_simulation_voltage();

      if (_voltages[it->first] != voltage)
        same_potentials = false;
        //accept_failure = false;

      _voltages[it->first] = voltage;

      if (voltage != 0.0)
        equilibrium = false;
    }
  }
  
  /* The same here: what about selfconsistency at 0 V ?? */
    //if (do_nothing)
    //  return;


  if (!equilibrium_done())
  {
    solve_equilibrium();

    // if we would repeat the equilibrium simulation, we can stop now
    if (equilibrium)
    {
      cout << "----------------------------------"
        << "--------------------------------->>" << endl;
      return;
    }
  }


  int coupling = get_options().coupling;
  
  if (equilibrium)
    get_options().coupling = POISSON;
  else if (same_potentials)
  {
    get_options().coupling = POISSON;
    try
    {
      do_newton();
    }
    catch (SolverException& e)
    {
      string msg = "DriftDiffusion: solve failed (" +
        string(e.what()) + ")";
      throw SolveFailedException(msg);
    }
    get_options().coupling = coupling;
  }

  /*
  static map<ElectricalContact*, double> voltages;
  ModelOptions& opts = SimulationInterface::get_options();
  bool quasi_equilibrium = false;
  if (opts.find_option("quasi_equilibrium"))
  {
    vector<string> qfpot(2, "");
    opts.get_option("quasi_equilibrium", qfpot);
    assert(qfpot.size() == 2);
    Boundary* boundary1 = get_environment().get_boundary(qfpot[0]);
    Boundary* boundary2 = get_environment().get_boundary(qfpot[1]);
    ElectricalContact* contact1;
    ElectricalContact* contact2;
    if ((boundary1 != NULL) && (boundary2 != NULL))
    {
      contact1 = dynamic_cast<ElectricalContact*>(
          boundary1->get_boundary_properties(get_id()));
      contact2 = dynamic_cast<ElectricalContact*>(
          boundary2->get_boundary_properties(get_id()));

      if ((contact1 != NULL) && (contact2 != NULL))
      {
        if ((contact1->get_simulation_voltage() != voltages[contact1]) ||
            (contact2->get_simulation_voltage() != voltages[contact2]))
        {
          voltages[contact1] = contact1->get_simulation_voltage();
          voltages[contact2] = contact2->get_simulation_voltage();
          set_electron_fermi_level(contact1->get_simulation_voltage());
          set_hole_fermi_level(contact2->get_simulation_voltage());
          get_options().coupling = POISSON;
          quasi_equilibrium = true;
          cerr << "solving quasi-equilibrium..." << endl;
        }
      }
    }
  }
  */

  if (do_local_scaling_)
    build_local_scaling();


  //set_dirichlet_bc();
  
  try
  {
    switch (_options.solver_method)
    {
      case GUMMEL:
        solve_gummel();
        break;
      default: // Newton method
        do_newton();
        break;
    }
  }
  catch (SolverException& e)
  {
    if (!accept_failure)
    {
      string msg = "DriftDiffusion: solve failed (" +
        string(e.what()) + ")";
      throw SolveFailedException(msg);
    }
  }

  get_options().coupling = coupling;

  // calculate the currents to print them on screen
  calculate_currents();

  ContactData::iterator it(_boundary_currents.begin());
  const ContactData::iterator end(_boundary_currents.end());

  cout << endl;
  int width = 20;
  {
    ostringstream os;
    os << "contact name:";
    os.width(width - os.tellp());
    os << "";
    os << "contact voltage:";
    os.width(2 * width - os.tellp());
    os << "";
    os << "contact current:";
    cout << os.str() << endl;
  }

  for ( ; it != end; ++it)
  {
    ostringstream os;
    os << setprecision(6);
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(it->first->get_boundary_properties(get_id()));
    os << it->first->get_name();
    os.width(width - os.tellp());
    os << "";
    os << cnt->get_simulation_voltage();
    os.width(2 * width - os.tellp());
    os << "";
    os << it->second * it->first->get_area_factor();
    cout << os.str() << endl;
  }
  cout << endl;
  cout << "------------------------------------------------------------------->>"
    << endl;

}




void
DriftDiffusion::do_equilibrium(void)
{

  // rebuild the system if needed
  //rebuild_equation_system();

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;
  
  //parse_options();


  // first we have to compute the scaling
  compute_scaling(get_options().scaling_type);

  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());
  
  SolverParameters& solver_params = get_options().solver_params;
  solver_params.nonlinear_max_iterations = 150;


  int coupling = get_options().coupling;
  get_options().coupling = POISSON;

  // backup the simulation voltages and set all to zero
  ContactData sim_voltages(_boundary_currents);
  ContactData::iterator it(sim_voltages.begin());
  const ContactData::iterator end(sim_voltages.end());
  for ( ; it != end; ++it)
  {
    const Boundary* bd = it->first;
    // It's save to static_cast because we know there has to be an
    // ElectricalContact object
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(bd->get_boundary_properties(get_id()));
    sim_voltages[bd] = cnt->get_simulation_voltage();
    cnt->set_simulation_voltage(0.0);
  }

  // make a rough guess
  guess_equilibrium();


  //if (do_local_scaling_)
  //  build_local_scaling();


  try
  {
    cerr << "Solving equilibrium" << endl;

    do_newton();
    
    cerr << "Equilibrium done" << endl;
  }
  catch (runtime_error& e)
  {
    cerr << "ATTENTION: Equilibrium did not converge: " << e.what() << endl;
    throw (e);
  }

  // set the contact voltages back to the desired values
  it = sim_voltages.begin();
  for ( ; it != end; ++it)
  {
    const Boundary* bd = it->first;
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(bd->get_boundary_properties(get_id()));
    cnt->set_simulation_voltage(sim_voltages[bd]);
  }
  
  // reset the coupling
  get_options().coupling = coupling;
}



void
DriftDiffusion::guess_equilibrium(void)
{

  // equation system needs to be active
  rebuild_equation_system();
  
  TiberNonlinearSystem& poisson =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const unsigned int u_var = poisson.variable_number("potential");
  const DofMap& dof_map_u = poisson.get_dof_map();
  vector<unsigned int> dof_indices_u;
  
  NumericVector<Number>& solution_u = poisson.get_solution_vector();
  solution_u.zero();

  MeshBase::const_element_iterator el =
                                  get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  get_mesh().active_elements_end();

  const double phi0 = get_scaling().get_potential_scaling();
  
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
      dynamic_cast<DriftDiffusionProperties*>(
          _device->get_material(elem->subdomain_id())->get_model(get_id()));

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
DriftDiffusion::parse_const_options(void)
{
  SolverParameters& solver_params = get_options().solver_params;

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

  string method =  opts.get_option("current_integration_method", "");
  if (method == "surfint")
    myopts.current_calculation = SURFINT;
  else
    myopts.current_calculation = RSTF;

  string scaling = opts.get_option("scaling", "");
  if (scaling == "demari")
    myopts.scaling_type = Scaling::DEMARI;
  else if (scaling == "none")
    myopts.scaling_type = Scaling::NONE;
  else
    myopts.scaling_type = Scaling::UNITS;

  do_local_scaling_ = opts.get_option("local_scaling", true);

  string discretization = opts.get_option("discretization", "");
  if (discretization == "fem")
    myopts.scheme = FEM;
  else if (discretization == "box")
    myopts.scheme = BOX;
  else if (discretization == "sg")
    myopts.scheme = SG;

  solver_params.nonlinear_solver =
    opts.get_option("nonlinear_solver", solver_params.nonlinear_solver);
  solver_params.linear_solver =
    opts.get_option("linear_solver", solver_params.linear_solver);

  string ksptype = opts.get_option("ksp_type", "");
  if (ksptype == "") {}
  else if (ksptype == "bcgsl")
    solver_params.ksp_type = BICGSTAB;
  else if (ksptype == "gmres")
    solver_params.ksp_type = GMRES;
  else if (ksptype == "bcgs")
    solver_params.ksp_type = BICG;
  else if (ksptype == "cg")
    solver_params.ksp_type = CG;
  else if (ksptype == "richardson")
    solver_params.ksp_type = RICHARDSON;



  string pc = opts.get_option("pc_type", "");
  if (pc == "") {}
  else if (pc == "ilu")
    solver_params.pc_type = ILU_PRECOND;
  else if (pc == "composite")
    solver_params.pc_type = USER_PRECOND;
  else if (pc == "jacobi")
    solver_params.pc_type = JACOBI_PRECOND;
  else if (pc == "lu")
    solver_params.pc_type = LU_PRECOND;
  else if (pc == "cholesky")
    solver_params.pc_type = CHOLESKY_PRECOND;
  else if (pc == "eisenstat")
    solver_params.pc_type = EISENSTAT_PRECOND;

}





void
DriftDiffusion::parse_options(void)
{
  PerfLog perf_log("DriftDiffusion parse_options()", false);
  perf_log.start_event("parse");
  
  SolverParameters& solver_params = get_options().solver_params;

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", 5));

  string coupling = opts.get_option("coupling", "");
  if (coupling == "full")
    myopts.coupling = FULLYCOUPLED;
  else if (coupling == "poisson")
    myopts.coupling = POISSON;
  else if (coupling == "electrons")
    myopts.coupling = ECURRENT | POISSON;
  else if (coupling == "holes")
    myopts.coupling = HCURRENT | POISSON;
  else if (coupling == "current")
    myopts.coupling = CURRENTS;
          

  myopts.mesh_refinement = opts.get_option("mesh_refinement",
      myopts.mesh_refinement);

  myopts.exact_newton = opts.get_option("exact_newton", myopts.exact_newton);

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
  solver_params.linear_max_iterations = opts.get_option("lin_max_it", 500);


  
  string lstype = opts.get_option("ls_type", "");
  if (lstype == "") {}
  else if (lstype == "none")
    solver_params.ls_type = 1;
  else if (lstype == "cubic")
    solver_params.ls_type = 3;
  else if (lstype == "quadratic")
    solver_params.ls_type = 2;

  solver_params.ls_maxstep = opts.get_option("ls_max_step",
      solver_params.ls_maxstep);

  perf_log.stop_event("parse");
}





void
DriftDiffusion::rebuild_equation_system(void)
{
  if (!_rebuild_eq_system) return;


  EquationSystems& equation_systems = get_equation_systems();

  SolverParameters& solver_params = get_options().solver_params;

  // the coupled DD system
  TiberNonlinearSystem& system =
    TiberNonlinearSystem::create_nonlinear_system(equation_systems,
        get_equation_system_name(), solver_params.nonlinear_solver,
        solver_params.linear_solver);

  
  const unsigned int dim = get_mesh().mesh_dimension();

  // TODO gives some problem in rev. 398
  // in 1D bcgs seems to work better than bcgsl
  if (dim == 1)
    if (solver_params.ksp_type == BICGSTAB)
      solver_params.ksp_type = BICG;

  system.set_linear_solver_type(solver_params.ksp_type, solver_params.pc_type);

  system.attach_assembly_routine(assemble_system);


  system.add_variable("potential", libMeshEnums::FIRST);
  system.add_variable("fermi_e", libMeshEnums::FIRST);
  system.add_variable("fermi_h", libMeshEnums::FIRST);

  // the linear solver is used to get good guesses for the electro-chemical
  // potentials
  _linear_solver = TiberLinearSolver::create("petsc");


  // finally initialize the newly created system
  system.init();

  _rebuild_eq_system = false;

}




void
DriftDiffusion::do_init(void)
{

  _device = &get_environment().get_device();

  find_dirichlet_nodes();
  find_dielectric_boundary_nodes();

  parse_const_options();
  
  rebuild_equation_system();

  // prepare the _boundary_currents
  // we will rely on the fact that it contains an entry for every boundary
  // later on !!!
  SimulationEnvironment::BoundaryIterator
    it(get_environment().boundaries_begin());
  const SimulationEnvironment::BoundaryIterator
    end(get_environment().boundaries_end());
  for ( ; it != end; ++it)
  {
    BoundaryProperties* bd = it->second->get_boundary_properties(get_id());
    if (bd != NULL)
    {
      _boundary_currents[it->second] = 0.0;
      _voltages[it->second] = 0.0;
    }
  }
}




NumericVector<double>&
DriftDiffusion::get_solution_vector(void)
{
  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());
  
  return system.get_solution_vector();
}




void
DriftDiffusion::do_newton(void) throw (SolverException)
{

  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());

  SolverParameters& solver_params = get_options().solver_params;
  
  system.set_linear_solver_params(solver_params.linear_tolerance,
      solver_params.linear_abs_tolerance,
      solver_params.linear_max_iterations);

  system.set_nonlinear_solver_params(solver_params.nonlinear_tolerance,
      solver_params.nonlinear_abs_tolerance,
      solver_params.nonlinear_step_tolerance,
      solver_params.ls_maxstep,
      solver_params.nonlinear_max_iterations);

  system.solve();
}







double
DriftDiffusion::do_gummel_iterations(int max_it) throw (SolverException)
{
 
  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  NumericVector<Number>& u = system.get_solution_vector();
  NumericVector<Number>& du = *(system.solution);
  
  SolverParameters& solver_params = get_options().solver_params;
  int nonlin_max_it = solver_params.nonlinear_max_iterations;
  int coupling = get_options().coupling;

  try
  {
    for (int i = 0; i < max_it; i++)
    {
      get_options().coupling = ECURRENT;
      do_newton();
      double norm_du = max(norm_du, du.linfty_norm());

      get_options().coupling = HCURRENT;
      do_newton();
      norm_du = max(norm_du, du.linfty_norm());

      get_options().coupling = POISSON;
      do_newton();
      norm_du = du.linfty_norm();

      double norm_res = system.rhs->l2_norm();
    cerr << "  Gummel it " << i << ", |du| = " << norm_du << ", |r| = " << norm_res << endl;

    }
  }
  catch (SolverException& err)
  {
    get_options().coupling = coupling;
    throw(err);
  }
  get_options().coupling = coupling;
 
  return 0;
}




void
DriftDiffusion::solve_gummel(void)
{
  do_gummel_iterations(get_options().max_gummel_iterations);
}




void
DriftDiffusion::get_band_edges(const Elem* elem, vector<double>& band_edges)
{
  band_edges.resize(2);

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
  }

  if (el != NULL) // we found it!
    get_bands_secure(el, band_edges);
  else
  {
    // no parent, so check for children
    vector<const Elem*> tree;
    elem->family_tree(tree, false);

    set<const Elem*> elem_list;
    unsigned int len = tree.size();
    for (unsigned int i = 0; i < len; i++)
    {
      const Elem* elem_i = tree[i];
      if (env.contains_element(elem_i))
        elem_list.insert(elem_i);
    }
    set<const Elem*>::iterator el_it = elem_list.begin();
    set<const Elem*>::iterator el_end = elem_list.end();
    for ( ; el_it != el_end; ++el_it)
    {
      el = *el_it;
      if (el->contains_point(elem->centroid()))
      {
        get_bands_secure(el, band_edges);
        // we have found it, so get out of the for loop
        break;
      }
    }
  }
}




void
DriftDiffusion::get_bands_secure(const Elem* elem, vector<double>& band_edges)
{

  assert(_device != NULL);

  DriftDiffusionProperties* sc =
    dynamic_cast<DriftDiffusionProperties*>(
        _device->get_material(elem->subdomain_id())->get_model(get_id()));
  assert(sc != NULL); 

  sc->reinit(elem);

  band_edges[0] = sc->get_conduction_band_edge();
  band_edges[1] = sc->get_valence_band_edge();

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

  SimulationEnvironment& env = get_environment();
  
  // check if elem is an active element of the simulation
  if (!env.contains_element(elem))
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();

    while ((parent != NULL) && (!env.contains_element(parent)))
      parent = parent->parent();

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
      const Elem* elem_i = tree[i];
      if (env.contains_element(elem_i))
        elem_list.insert(elem_i);
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
      if (el_it == el_end)
        solution[i] = 0;
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
  else
    solution = 0;
}


ID
DriftDiffusion::convert_variable_name_to_id(const string& variable_name) const
{
  ID id = INVALID_ID;

  // for an empty string we return immediately
  if (variable_name == "") return id;

  switch (variable_name[0])
  {
    case 'E':
      if (variable_name == "ElPotential")
        id = ELPOTENTIAL;
      else if (variable_name == "Ev")
        id = VBANDEDGE;
      else if (variable_name == "Ec")
        id = CBANDEDGE;
      else if (variable_name == "Eg")
        id = BANDGAP;
      else if (variable_name == "Ev0")
        id = VBANDEDGEINTR;
      else if (variable_name == "Ec0")
        id = CBANDEDGEINTR;
      break;

    case 'e':
      if (variable_name == "eDensity")
        id = EDENSITY;
      else if (variable_name == "eMob")
        id = EMOBILITY;
      else if (variable_name == "eCond")
        id = SIGMAE;
      break;
      
    case 'h':
      if (variable_name == "hDensity")
        id = HDENSITY;
      else if (variable_name == "hMob")
        id = HMOBILITY;
      else if (variable_name == "hCond")
        id = SIGMAH;
      break;
      
    case 'J':
      if (variable_name.size() == 1)
        id = J;
      else
      {
        switch (variable_name[1])
        {
          case 'n':
            if (variable_name.size() == 2)
              id = JN;
            else if (variable_name == "Jn_x")
              id = JNX;
            else if (variable_name == "Jn_y")
              id = JNY;
            if (variable_name == "Jn_z")
              id = JNZ;
            break;
          case 'p':
            if (variable_name.size() == 2)
              id = JP;
            else if (variable_name == "Jp_x")
              id = JPX;
            else if (variable_name == "Jp_y")
              id = JPY;
            if (variable_name == "Jp_z")
              id = JPZ;
            break;
        }
      }
      break;
      
    case 'Q':
      if (variable_name == "QFermi_e")
        id = QFERMIE;
      else if (variable_name == "QFermi_h")
        id = QFERMIH;
      break;

    case 'r':
      {
        vector<string> rec;
        Utils::tokenize(variable_name, rec);
        if (rec.size() > 1)
        {
          if (rec[0] == "recombination")
          {
            if (rec[1] == "total")
              id = RECOMB;
            else
            {
              ID rec_id = PhysicalModelInterface::get_id_from_name<
                RecombinationModelInterface>(rec[1]);
              id = MODELS + rec_id;
            }
          }
        }
      }
      

    default:
      break;
  }

  return id;
}

      


void
DriftDiffusion::get_solution_secure(const Elem* elem,
    const set<ID>& ids, vector<map<ID, double> >& values)
{
  unsigned int np = elem->n_nodes();

  vector<Point> points(np);
  for (int i = 0; i < np; i++)
    points[i] = elem->local_node(elem->type(), i);

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const Device& device = *(_device);

  const NumericVector<Number>& solution = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  //FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);


  ID subdomain = elem->subdomain_id();

  DriftDiffusionProperties* sc =
    dynamic_cast<DriftDiffusionProperties*>(
        device.get_material(subdomain)->get_model(get_id()));
  assert(sc != NULL); 

  sc->reinit(elem);

  fe->reinit(elem, &points);
    
  //Get the thermoelectric power
  double Pn =  get_electrons_thermoelectric_power(elem);
  double Pp =  get_holes_thermoelectric_power(elem);
 
  vector<double> T_nodes = sc->get_temperature_node();

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  for (unsigned int n = 0; n < np; n++)
  {
    double en_x = 0.0, ep_x = 0.0;
    double en_y = 0.0, ep_y = 0.0;
    double en_z = 0.0, ep_z = 0.0;
    double dT_x = 0.0;
    double dT_y = 0.0;
    double dT_z = 0.0;
    double u  = phi0 * solution(dof_indices_u[n]);
    double en = phi0 * solution(dof_indices_en[n]);
    double ep = phi0 * solution(dof_indices_ep[n]);
    RealGradient e_field(0);
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      en_x  += dphi[i][n](0) * solution(dof_indices_en[i]);
      en_y  += dphi[i][n](1) * solution(dof_indices_en[i]);
      en_z  += dphi[i][n](2) * solution(dof_indices_en[i]);
      
      ep_x  += dphi[i][n](0) * solution(dof_indices_ep[i]);
      ep_y  += dphi[i][n](1) * solution(dof_indices_ep[i]);
      ep_z  += dphi[i][n](2) * solution(dof_indices_ep[i]);
      
      dT_x  += dphi[i][n](0) * T_nodes[i];
      dT_y  += dphi[i][n](1) * T_nodes[i];
      dT_z  += dphi[i][n](2) * T_nodes[i];

      e_field += dphi[i][n] * solution(dof_indices_u[i]);
    }

    // scale the potential back
    e_field *= -phi0;
    en_x *= phi0;
    en_y *= phi0;
    en_z *= phi0;
    ep_x *= phi0;
    ep_y *= phi0;
    ep_z *= phi0;

    sc->set_coordinates(points[n]);
    // all are at lattice temperature

    sc->set_potentials(u, en, ep);

    sc->set_electric_field(e_field);

    sc->calculate_densities();
    sc->calculate_mobilities();

    // we put the minus here for convenience
    double sigma_e = -Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = -Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();

    double jnx = sigma_e * (en_x + Pn * dT_x);
    double jny = sigma_e * (en_y + Pn * dT_y);
    double jnz = sigma_e * (en_z + Pn * dT_z);
    double jpx = sigma_h * (ep_x + Pp * dT_x);
    double jpy = sigma_h * (ep_y + Pp * dT_y);
    double jpz = sigma_h * (ep_z + Pp * dT_z);


    if (ids.count(ELPOTENTIAL))
      values[n][ELPOTENTIAL] = u;

    if (ids.count(QFERMIE))
      values[n][QFERMIE] = en;

    if (ids.count(QFERMIH))
      values[n][QFERMIH] = ep;

    if (ids.count(VBANDEDGE))
      values[n][VBANDEDGE] = sc->get_valence_band_edge() - u;

    if (ids.count(CBANDEDGE))
      values[n][CBANDEDGE] = sc->get_conduction_band_edge() - u;

    if (ids.count(VBANDEDGEINTR))
      values[n][VBANDEDGEINTR] = sc->get_valence_band_edge();

    if (ids.count(CBANDEDGEINTR))
      values[n][CBANDEDGEINTR] = sc->get_conduction_band_edge();

    if (ids.count(EDENSITY))
      values[n][EDENSITY] = sc->get_electron_density();

    if (ids.count(HDENSITY))
      values[n][HDENSITY] = sc->get_hole_density();

    if (ids.count(BANDGAP))
      values[n][BANDGAP] =
        sc->get_conduction_band_edge() - sc->get_valence_band_edge();

    if (ids.count(EMOBILITY))
      values[n][EMOBILITY] = sc->get_electron_mobility();

    if (ids.count(HMOBILITY))
      values[n][HMOBILITY] = sc->get_hole_mobility();

    if (ids.count(SIGMAE))
      values[n][SIGMAE] = Constants::e *
        sc->get_electron_mobility() * sc->get_electron_density();

    if (ids.count(SIGMAH))
      values[n][SIGMAH] = Constants::e *
        sc->get_hole_mobility() * sc->get_hole_density();

    if (ids.count(J))
    {
      double jx = jnx + jpx;
      double jy = jny + jpy;
      double jz = jnz + jpz;
      values[n][J] = sqrt(jx * jx + jy * jy + jz * jz);
    }

    if (ids.count(JN))
      values[n][JN] = sqrt(jnx * jnx + jny * jny + jnz * jnz);

    if (ids.count(JP))
      values[n][JP] = sqrt(jpx * jpx + jpy * jpy + jpz * jpz);

    if (ids.count(JNX))
      values[n][JNX] = jnx;

    if (ids.count(JNY))
      values[n][JNY] = jny;

    if (ids.count(JNZ))
      values[n][JNZ] = jnz;

    if (ids.count(JPX))
      values[n][JPX] = jpx;

    if (ids.count(JPY))
      values[n][JPY] = jpy;

    if (ids.count(JPZ))
      values[n][JPZ] = jpz;

    set<ID>::iterator first(ids.begin());
    set<ID>::iterator it(ids.end());
    --it;
    while (*it > MODELS)
    {
      values[n][*it] = sc->get_net_recombination_rate(*it - MODELS);
      if (it == first)
        break;

      --it;
    }
  }
}

      


void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{
  unsigned int np = p.size();

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const Device& device = *(_device);

  const NumericVector<Number>& solution = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);


  ID subdomain = elem->subdomain_id();

  DriftDiffusionProperties* sc =
    dynamic_cast<DriftDiffusionProperties*>(
        device.get_material(subdomain)->get_model(get_id()));
  assert(sc != NULL); 

  sc->reinit(elem);

  fe->reinit(elem, &points);
    
  //Get the thermoelectric power
  double Pn =  get_electrons_thermoelectric_power(elem);
  double Pp =  get_holes_thermoelectric_power(elem);
 
  vector<double> T_nodes = sc->get_temperature_node();

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  for (unsigned int n = 0; n < np; n++)
  {
    double en_x = 0.0, ep_x = 0.0;
    double en_y = 0.0, ep_y = 0.0;
    double en_z = 0.0, ep_z = 0.0;
    double dT_x = 0.0;
    double dT_y = 0.0;
    double dT_z = 0.0;
    double u = 0;
    double en = 0.0;
    double ep = 0.0;
    RealGradient e_field(0);
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u  += phi[i][n] * solution(dof_indices_u[i]);
      en += phi[i][n] * solution(dof_indices_en[i]);
      ep += phi[i][n] * solution(dof_indices_ep[i]);

      en_x  += dphi[i][n](0) * solution(dof_indices_en[i]);
      en_y  += dphi[i][n](1) * solution(dof_indices_en[i]);
      en_z  += dphi[i][n](2) * solution(dof_indices_en[i]);
      
      ep_x  += dphi[i][n](0) * solution(dof_indices_ep[i]);
      ep_y  += dphi[i][n](1) * solution(dof_indices_ep[i]);
      ep_z  += dphi[i][n](2) * solution(dof_indices_ep[i]);
      
      dT_x  += dphi[i][n](0) * T_nodes[i];
      dT_y  += dphi[i][n](1) * T_nodes[i];
      dT_z  += dphi[i][n](2) * T_nodes[i];

      e_field += dphi[i][n] * solution(dof_indices_u[i]);
    }

    // scale the potential back
    u  *= phi0;
    en  *= phi0;
    ep  *= phi0;
    e_field *= -phi0;
    en_x *= phi0;
    en_y *= phi0;
    en_z *= phi0;
    ep_x *= phi0;
    ep_y *= phi0;
    ep_z *= phi0;

    sc->set_coordinates(points[n]);

    sc->set_potentials(u, en, ep);

    sc->set_electric_field(e_field);

    sc->calculate_densities();
    sc->calculate_mobilities();

    // we put the minus here for convenience
    double sigma_e = -Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = -Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();

    double jnx = sigma_e * (en_x + Pn * dT_x);
    double jny = sigma_e * (en_y + Pn * dT_y);
    double jnz = sigma_e * (en_z + Pn * dT_z);
    double jpx = sigma_h * (ep_x + Pp * dT_x);
    double jpy = sigma_h * (ep_y + Pp * dT_y);
    double jpz = sigma_h * (ep_z + Pp * dT_z);


    if (ids.count(ELPOTENTIAL))
      values[n][ELPOTENTIAL] = u;

    if (ids.count(QFERMIE))
      values[n][QFERMIE] = en;

    if (ids.count(QFERMIH))
      values[n][QFERMIH] = ep;

    if (ids.count(VBANDEDGE))
      values[n][VBANDEDGE] = sc->get_valence_band_edge() - u;

    if (ids.count(CBANDEDGE))
      values[n][CBANDEDGE] = sc->get_conduction_band_edge() - u;

    if (ids.count(VBANDEDGEINTR))
      values[n][VBANDEDGEINTR] = sc->get_valence_band_edge();

    if (ids.count(CBANDEDGEINTR))
      values[n][CBANDEDGEINTR] = sc->get_conduction_band_edge();

    if (ids.count(EDENSITY))
      values[n][EDENSITY] = sc->get_electron_density();

    if (ids.count(HDENSITY))
      values[n][HDENSITY] = sc->get_hole_density();

    if (ids.count(BANDGAP))
      values[n][BANDGAP] =
        sc->get_conduction_band_edge() - sc->get_valence_band_edge();

    if (ids.count(EMOBILITY))
      values[n][EMOBILITY] = sc->get_electron_mobility();

    if (ids.count(HMOBILITY))
      values[n][HMOBILITY] = sc->get_hole_mobility();

    if (ids.count(J))
    {
      double jx = jnx + jpx;
      double jy = jny + jpy;
      double jz = jnz + jpz;
      values[n][J] = sqrt(jx * jx + jy * jy + jz * jz);
    }

    if (ids.count(JN))
      values[n][JN] = sqrt(jnx * jnx + jny * jny + jnz * jnz);

    if (ids.count(JP))
      values[n][JP] = sqrt(jpx * jpx + jpy * jpy + jpz * jpz);

    if (ids.count(JNX))
      values[n][JNX] = jnx;

    if (ids.count(JNY))
      values[n][JNY] = jny;

    if (ids.count(JNZ))
      values[n][JNZ] = jnz;

    if (ids.count(JPX))
      values[n][JPX] = jpx;

    if (ids.count(JPY))
      values[n][JPY] = jpy;

    if (ids.count(JPZ))
      values[n][JPZ] = jpz;

    set<ID>::iterator first(ids.begin());
    set<ID>::iterator it(ids.end());
    --it;
    while (*it > MODELS)
    {
      values[n][*it] = sc->get_net_recombination_rate(*it - MODELS);
      if (it == first)
        break;

      --it;
    }

  }
}




template <>
void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<double>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& ddsol = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

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





template <>
void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<DriftDiffusion::Solution>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& ddsol = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

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


template <>
void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<DriftDiffusion::Currents>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& ddsol = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);

  ID subdomain = elem->subdomain_id();

  DriftDiffusionProperties* sc =
    dynamic_cast<DriftDiffusionProperties*>(
        _device->get_material(subdomain)->get_model(get_id()));

  assert(sc != NULL); 

  sc->reinit(elem);


  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  const double phi0 = get_scaling().get_potential_scaling();


  for (unsigned int n = 0; n < np; n++)
  {
    double u = 0;
    double en = 0;
    double ep = 0;
    double en_x = 0.0, ep_x = 0.0;
    double en_y = 0.0, ep_y = 0.0;
    double en_z = 0.0, ep_z = 0.0;
    RealGradient e_field(0);
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u  += phi[i][n] * ddsol(dof_indices_u[i]);
      en += phi[i][n] * ddsol(dof_indices_en[i]);
      ep += phi[i][n] * ddsol(dof_indices_ep[i]);

      en_x  += dphi[i][n](0) * ddsol(dof_indices_en[i]);
      en_y  += dphi[i][n](1) * ddsol(dof_indices_en[i]);
      en_z  += dphi[i][n](2) * ddsol(dof_indices_en[i]);
      
      ep_x  += dphi[i][n](0) * ddsol(dof_indices_ep[i]);
      ep_y  += dphi[i][n](1) * ddsol(dof_indices_ep[i]);
      ep_z  += dphi[i][n](2) * ddsol(dof_indices_ep[i]);

      e_field += dphi[i][0] * ddsol(dof_indices_u[i]);
    }

    // scale the potential back
    e_field *= -phi0;
    en_x *= phi0;
    en_y *= phi0;
    en_z *= phi0;
    ep_x *= phi0;
    ep_y *= phi0;
    ep_z *= phi0;

    // prepare for calculating local properties
    sc->set_coordinates(points[n]);


    sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

    sc->set_electric_field(e_field);

    sc->calculate_densities();
    sc->calculate_mobilities();

    // we put the minus here for convenience
    double sigma_e = -Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = -Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();


    solution[n]._jn_x = sigma_e * en_x;
    solution[n]._jn_y = sigma_e * en_y;
    solution[n]._jn_z = sigma_e * en_z;

    solution[n]._jp_x = sigma_h * ep_x;
    solution[n]._jp_y = sigma_h * ep_y;
    solution[n]._jp_z = sigma_h * ep_z;
  }
}


template <>
void
DriftDiffusion::get_solution_secure(const Elem* elem, const vector<Point>& p,
    vector<DriftDiffusion::EField>& solution)
{
  unsigned int np = p.size();
  solution.resize(np);
  if (np == 0) return;

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& ddsol = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  vector<Point> points(np);
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
  //for (unsigned int n = 0; n < np; n++)
  //  points[n] = FEInterface::inverse_map(dim, fe_type, elem, p[n], 1e-6);

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_u, u_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  const double phi0 = get_scaling().get_potential_scaling();


  for (unsigned int n = 0; n < np; n++)
  {
    RealGradient& efield = solution[n]._efield;

    efield = 0;
    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
      efield += dphi[i][0] * ddsol(dof_indices_u[i]);

    // scale the potential back
    efield *= -phi0;
  }
}







void
DriftDiffusion::get_solution(const Elem* elem,
    std::vector<DriftDiffusion::Solution>& solution)
{
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const DofMap& dof_map = system->get_dof_map();
  const NumericVector<Number>& ddsol = system->get_solution_vector();

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
    solution[n].potential = phi0 * ddsol(dof_indices_u[n]);
    solution[n].fermi_e = phi0 * ddsol(dof_indices_en[n]);
    solution[n].fermi_h = phi0 * ddsol(dof_indices_ep[n]);
  }

}






void
DriftDiffusion::calculate_currents_rstf(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
  // reset currents
  ContactData::iterator it =
    _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second = 0.0;

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const Mesh& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();


  const double phi0 = get_scaling().get_potential_scaling();

  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  FEType fe_type = system->variable_type(u_var);

  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  QGauss qrule(dim, get_options().integration_order);
  fe->attach_quadrature_rule(&qrule);

  
  // Jacobian * quadrature weight at each integration point.   
  const vector<Real>& JxW = fe->get_JxW();
  
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;


  // will contain the node ids if an element has boundary nodes
  vector<Boundary*> node_ids;

  MeshBase::const_element_iterator el =
                                  mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();
        
    bool has_node = false;
    node_ids.resize(elem->n_nodes());
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      Boundary* bd = env.get_boundary(elem->get_node(n));
      node_ids[n] = bd;
      if (bd != NULL)
        has_node = true;
    }

    // if the element has no node on a boundary,
    // we can go to the next element
    if (!has_node)
      continue;
    
    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);

    
    // in a dielectric we have no current
    if (sc->is_dielectric())
      continue;


    fe->reinit(elem);

    sc->reinit(elem);
    
    //Get the thermoelectric power------------
    double Pn =  get_electrons_thermoelectric_power(elem) / phi0;
    double Pp =  get_holes_thermoelectric_power(elem) / phi0;

    //Get the temperature given the element
    vector<double> T_nodes =  sc->get_temperature_node();

        
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      unsigned int n_dofs = dof_indices_u.size();
      // get the solution values at the centroid
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      RealGradient dEfn(0);
      RealGradient dEfp(0);
      RealGradient e_field(0);
      RealGradient dT(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        ep += phi[i][qp] * solution(dof_indices_ep[i]);

        dEfn += dphi[i][qp] * solution(dof_indices_en[i]);
        dEfp += dphi[i][qp] * solution(dof_indices_ep[i]);

        dT += dphi[i][qp] * T_nodes[i];

        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }

      // prepare for calculating local properties
      sc->set_coordinates(elem->centroid());


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

      sc->set_electric_field(phi0 * e_field);

      sc->calculate_densities();
      sc->calculate_mobilities();

      // we put the minus here for convenience
      double sigma_e = -Constants::e * sc->get_electron_density() *
        sc->get_electron_mobility();
      double sigma_h = -Constants::e * sc->get_hole_density() *
        sc->get_hole_mobility();

      RealGradient j(JxW[qp] * phi0 *
          (sigma_e * (dEfn + Pn * dT) + sigma_h * (dEfp + Pp * dT))); 

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        Boundary* boundary = node_ids[n];
        if (boundary != NULL)
          _boundary_currents[boundary] += j * dphi[n][qp];

      }
    } // end loop over quadrature points
  } // end loop over elements

}




void
DriftDiffusion::calculate_currents_surfint(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;
  
  // reset currents
  ContactData::iterator it =
    _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second = 0.0;

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const Mesh& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = get_scaling().get_potential_scaling();

  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type));
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

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(top_parent, s);
      
      if (env.is_boundary(side))
      {

        Boundary* boundary = env.get_boundary(side);
        if (boundary == NULL)
          continue;
        
        sc->reinit(elem);
        
	//Get the thermoelectric power------------
        double Pn = get_electrons_thermoelectric_power(elem) / phi0;
        double Pp = get_holes_thermoelectric_power(elem) / phi0;
       
        //Get the temperature given the element
        vector<double> T_nodes = sc->get_temperature_node();

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
            RealGradient e_field(0);
            Real dT = 0.0;
            for (unsigned int i = 0; i < phi_size; i++)
            {
              u  += phi[i][qp] * solution(dof_indices_u[i]);
              en += phi[i][qp] * solution(dof_indices_en[i]);
              ep += phi[i][qp] * solution(dof_indices_ep[i]);
              
              double tmp = dphi[i][qp] * face_normals[qp];
              dEfn += tmp * solution(dof_indices_en[i]);
              dEfp += tmp * solution(dof_indices_ep[i]);

              e_field += dphi[i][qp] * solution(dof_indices_u[i]);

              dT += tmp * T_nodes[i];
            }

            // prepare for calculating local properties
            sc->set_coordinates(q_point[qp]);


            sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

            sc->set_electric_field(phi0 * e_field);

            sc->calculate_densities();

            sc->calculate_mobilities();


            Real cond_e = Constants::e * sc->get_electron_mobility() *
              sc->get_electron_density();
            Real cond_h = Constants::e * sc->get_hole_mobility() *
              sc->get_hole_density();

            current += -JxW[qp] * phi0 *
              (cond_e * (dEfn + Pn * dT) + cond_h * (dEfp + Pp * dT));
          } // end loop over quadrature points

          _boundary_currents[boundary] += current;
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
          RealGradient e_field(0);
          Real dT = 0.0;
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            dEfn += dphi[n][0](0) * solution(dof_indices_en[n]);
            dEfp += dphi[n][0](0) * solution(dof_indices_ep[n]);
            e_field(0) += dphi[n][0](0) * solution(dof_indices_u[n]);

            dT += dphi[n][0](0) * T_nodes[n];
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
            dT = -dT;
          }

          // prepare for calculating local properties
          sc->set_coordinates(elem->point(s));


          sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

          sc->set_electric_field(phi0 * e_field);

          sc->calculate_densities();

          sc->calculate_mobilities();


          Real cond_e = Constants::e * sc->get_electron_mobility() *
            sc->get_electron_density();
          Real cond_h = Constants::e * sc->get_hole_mobility() *
            sc->get_hole_density();

          _boundary_currents[boundary] = -phi0 *
            (cond_e * (dEfn + Pn * dT) + cond_h * (dEfp + Pp * dT));
        }
      }
    } // end loop over elem sides
  } // end loop over elements

}





void
DriftDiffusion::build_local_scaling(void)
{
  
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());
  
  const NumericVector<Number>& solution = system->get_solution_vector();

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
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;


  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  FEType fe_type = system->variable_type(u_var);
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




  local_scaling_.clear();
  {
    MeshBase::const_element_iterator it =
      mesh.active_elements_begin();
    const MeshBase::const_element_iterator end =
      mesh.active_elements_end(); 

    for ( ; it != end; ++it)
      for (unsigned int n = 0; n < (*it)->n_nodes(); n++)
        local_scaling_[(*it)->get_node(n)] = vector<double>(3, 0.0);
  }


  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end(); 

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    unsigned int n_dofs     = dof_indices_u.size();

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));
    assert(sc != NULL); 

    sc->reinit(elem);

    fe->reinit(elem);

    assert(elem->n_nodes() == dof_indices_u.size());

    RealGradient field(0.0);
    for (unsigned int i = 0; i < dof_indices_u.size(); i++)
      field += dphi[i][0] * solution(dof_indices_u[i]);
    field *= -phi0;


    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      RealGradient e_field(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        ep += phi[i][qp] * solution(dof_indices_ep[i]);
        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }



      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
      sc->set_electric_field(phi0 * e_field);

      sc->calculate_densities();
      sc->calculate_mobilities();

      double epsilon = sc->get_relative_permittivity();
      double l2_eps = JxW[qp] * l2 * epsilon;

      double sigma_e = JxW[qp] * //sc->get_electron_density();
        sc->get_electron_mobility() * sc->get_electron_density();
      double sigma_h = JxW[qp] * //sc->get_hole_density();
        sc->get_hole_mobility() * sc->get_hole_density();

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //local_scaling_[elem->get_node(i)][0] += sigma_e * phi[i][qp];
        //local_scaling_[elem->get_node(i)][1] += sigma_h * phi[i][qp];
        local_scaling_[elem->get_node(i)][0] +=
          sigma_e * (dphi[i][qp] * dphi[i][qp]);
        local_scaling_[elem->get_node(i)][1] +=
          sigma_h * (dphi[i][qp] * dphi[i][qp]);
        local_scaling_[elem->get_node(i)][2] +=
          l2_eps * (dphi[i][qp] * dphi[i][qp]);
      }


    } // end loop over quadrature points
  } // end loop over elements
}






// implementation taken from libmesh equation_systems.C
void
DriftDiffusion::build_nodal_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{
  
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const Device& device = *(_device);
  const Mesh& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();
  
  legend.resize(variables.size());

  // for each possible variable we set the vector index
  // -1 means, the variable should not be plotted
  unsigned int n_vars = 0;
  const set<string>::const_iterator varend(variables.end());

  int Ec = -1;
  if (variables.find("Ec") != varend)
  {
    Ec = n_vars;
    legend[n_vars] = "Ec";
    n_vars++;
  }

  int Ev = -1;
  if (variables.find("Ev") != varend)
  {
    Ev = n_vars;
    legend[n_vars] = "Ev";
    n_vars++;
  }

  int Ec0 = -1;
  if (variables.find("Ec0") != varend)
  {
    Ec0 = n_vars;
    legend[n_vars] = "Ec0";
    n_vars++;
  }

  int Ev0 = -1;
  if (variables.find("Ev0") != varend)
  {
    Ev0 = n_vars;
    legend[n_vars] = "Ev0";
    n_vars++;
  }


  int Efn = -1;
  if (variables.find("QFermi_e") != varend)
  {
    Efn = n_vars;
    legend[n_vars] = "QFermi_e";
    n_vars++;
  }

  int Efp = -1;
  if (variables.find("QFermi_h") != varend)
  {
    Efp = n_vars;
    legend[n_vars] = "QFermi_h";
    n_vars++;
  }

  int phi = -1;
  if (variables.find("ElPotential") != varend)
  {
    phi = n_vars;
    legend[n_vars] = "electric_potential";
    n_vars++;
  }

  int Eg = -1;
  if (variables.find("Eg") != varend)
  {
    Eg = n_vars;
    legend[n_vars] = "Eg";
    n_vars++;
  }

  int edens = -1;
  if (variables.find("eDensity") != varend)
  {
    edens = n_vars;
    legend[n_vars] = "electron_density";
    n_vars++;
  }

  int hdens = -1;
  if (variables.find("hDensity") != varend)
  {
    hdens = n_vars;
    legend[n_vars] = "hole_density";
    n_vars++;
  }

  int Nd = -1;
  if (variables.find("Nd") != varend)
  {
    Nd = n_vars;
    legend[n_vars] = "ionized_donors";
    n_vars++;
  }

  int Na = -1;
  if (variables.find("Na") != varend)
  {
    Na = n_vars;
    legend[n_vars] = "ionized_acceptors";
    n_vars++;
  }

  int rho = -1;
  if (variables.find("charge_density") != varend)
  {
    rho = n_vars;
    legend[n_vars] = "total_charge_densitity";
    n_vars++;
  }

  int rec = -1;
  int num_rec = 0;
  map<ID, string> rec_model_ids;
  if (variables.find("NetRecombination") != varend)
  {

    // look for all recombination models
    {
      MeshBase::const_element_iterator it =
        mesh.active_local_elements_begin();

      assert(it != mesh.active_local_elements_end());

      DriftDiffusionProperties* sc = NULL;
      vector<ID> ids;

      const Elem* elem = *it;

      ID subdomain = elem->subdomain_id();

      sc = dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

      int n = sc->get_net_recombination_rate_IDs(ids);

      for (int i = 0; i < n; i++)
        rec_model_ids[ids[i]] =
          sc->get_recombination_model(ids[i])->get_default_name();
    }
    num_rec = rec_model_ids.size();
    legend.resize(variables.size() + num_rec);
    if (num_rec != 0)
    {
      rec = n_vars;

      map<ID, string>::iterator it = rec_model_ids.begin();
      map<ID, string>::iterator end = rec_model_ids.end();
      for ( ; it != end; ++it)
      {
        legend[n_vars] = it->second;
        n_vars++;
      }

      if (num_rec > 1)
      {
        // plot also the total net rate
        legend[n_vars] = "total_net_recombination";
        n_vars++;
      }
    }
  }


  int mun = -1;
  if (variables.find("eMob") != varend)
  {
    mun = n_vars;
    legend[n_vars] = "electron_mobility";
    n_vars++;
  }

  int mup = -1;
  if (variables.find("hMob") != varend)
  {
    mup = n_vars;
    legend[n_vars] = "hole_mobility";
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

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

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
      dof_map.dof_indices(elem, dof_indices_en, en_var);
      dof_map.dof_indices(elem, dof_indices_ep, ep_var);

      DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));
      assert(sc != NULL); 

      sc->reinit(elem);

      fe->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      RealGradient field(0.0);
      for (unsigned int i = 0; i < dof_indices_u.size(); i++)
        field += dphi[i][0] * solution(dof_indices_u[i]);
      field *= -phi0;


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        double u  = phi0 * solution(dof_indices_u[n]);
        double en = phi0 * solution(dof_indices_en[n]);
        double ep = phi0 * solution(dof_indices_ep[n]);


        // prepare for calculating local properties
        sc->set_coordinates(elem->point(n));


        sc->set_potentials(u, en, ep);
        sc->set_electric_field(field);

        sc->calculate_densities();
        sc->calculate_ionized_dopants();
        sc->calculate_mobilities();
        sc->calculate_net_recombination_rates();


        assert (node_conn[elem->node(n)] != 0);
        double conn = static_cast<double>(node_conn[elem->node(n)]);

        unsigned int id = n_vars * elem->node(n);

        if (edens != -1)
        {
          double nodal_val = sc->get_electron_density();
          if (sc->is_dielectric())
            nodal_val = 0.0;
          local[id + edens] += nodal_val / conn;
        }

        if (hdens != -1)
        {
          double nodal_val = sc->get_hole_density();
          if (sc->is_dielectric())
            nodal_val = 0.0;
          local[id + hdens] += nodal_val / conn;
        }

        if (Nd != -1)
        {
          double nodal_val = sc->get_ionized_donor_density();
          local[id + Nd] += nodal_val / conn;
        }

        if (Na != -1)
        {
          double nodal_val = sc->get_ionized_acceptor_density();
          local[id + Na] += nodal_val / conn;
        }

        if (rho != -1)
        {
          double nodal_val = sc->get_charge_density();
          if (sc->is_dielectric())
            nodal_val = 0.0;
          local[id + rho] += nodal_val / conn;
        }

        if (rec != -1)
        {
          // recombination models
          map<ID, string>::iterator it = rec_model_ids.begin();
          map<ID, string>::iterator end = rec_model_ids.end();
          int ctr = id + rec;
          double tot = 0.0;
          for ( ; it != end; ++it, ctr++)
          {
            double nodal_val = sc->get_net_recombination_rate(it->first);
            local[ctr] += nodal_val / conn;
            tot += nodal_val;
          }
          if (num_rec > 1)
            local[ctr] += tot / conn;
        }


        if (mun != -1)
        {
          double nodal_val = sc->get_electron_mobility();
          local[id + mun] += nodal_val / conn;
        }

        if (mup != -1)
        {
          double nodal_val = sc->get_hole_mobility();
          local[id + mup] += nodal_val / conn;
        }


        if (phi != -1)
          local[id + phi] += u / conn;

        if (Efn != -1)
          local[id + Efn] += -en / conn;

        if (Efp != -1)
          local[id + Efp] += -ep / conn;

        if (Ec != -1)
          local[id + Ec] += (sc->get_conduction_band_edge() - u) / conn;

        if (Ev != -1)
          local[id + Ev] += (sc->get_valence_band_edge() - u) / conn;

        if (Ec0 != -1)
          local[id + Ec0] += sc->get_conduction_band_edge() / conn;

        if (Ev0 != -1)
          local[id + Ev0] += sc->get_valence_band_edge() / conn;


        if (Eg != -1)
          local[id + Eg] +=
            (sc->get_conduction_band_edge() - sc->get_valence_band_edge()) / conn;


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
DriftDiffusion::build_elemental_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{

  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;
  
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());
  
  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const Device& device = *(_device);
  const Mesh& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_active_elem();
  
  legend.resize(variables.size());

  // for each possible variable we set the vector index
  // -1 means, the variable should not be plotted
  unsigned int n_vars = 0;
  const set<string>::const_iterator varend(variables.end());

  int EField = -1;
  if (variables.find("EField") != varend)
  {
    legend.resize(legend.size() + dim);
    EField = n_vars;
    switch (dim)
    {
      case 3:
        legend[EField + 2] = "E_z";
        n_vars++;
      case 2:
        legend[EField + 1] = "E_y";
        n_vars++;
        legend[EField + dim] = "modE";
        n_vars++;
      default:
        legend[EField] = "E_x";
        n_vars++;
    }
  }

  int Jn = -1;
  if (variables.find("eCurrent") != varend)
  {
    legend.resize(legend.size() + dim);
    Jn = n_vars;
    switch (dim)
    {
      case 3:
        legend[Jn + 2] = "Jn_z";
        n_vars++;
      case 2:
        legend[Jn + 1] = "Jn_y";
        n_vars++;
        legend[Jn + dim] = "modJn";
        n_vars++;
      default:
        legend[Jn] = "Jn_x";
        n_vars++;
    }
  }


  int Jp = -1;
  if (variables.find("hCurrent") != varend)
  {
    legend.resize(legend.size() + dim);
    Jp = n_vars;
    switch (dim)
    {
      case 3:
        legend[Jp + 2] = "Jp_z";
        n_vars++;
      case 2:
        legend[Jp + 1] = "Jp_y";
        n_vars++;
        legend[Jp + dim] = "modJp";
        n_vars++;
      default:
        legend[Jp] = "Jp_x";
        n_vars++;
    }
  }


  int J = -1;
  if (variables.find("Current") != varend)
  {
    legend.resize(legend.size() + dim);
    J = n_vars;
    switch (dim)
    {
      case 3:
        legend[J + 2] = "J_z";
        n_vars++;
      case 2:
        legend[J + 1] = "J_y";
        n_vars++;
        legend[J + dim] = "modJ";
        n_vars++;
      default:
        legend[J] = "J_x";
        n_vars++;
    }
  }

  int Polariz = -1;
  if ((variables.find("Polarization") != varend) ||
      (variables.find("polarization") != varend))
  {
    legend.resize(legend.size() + dim);
    Polariz = n_vars;
    switch (dim)
    {
      case 3:
        legend[Polariz + 2] = "P_z";
        n_vars++;
      case 2:
        legend[Polariz + 1] = "P_y";
        n_vars++;
        legend[Polariz + dim] = "modP";
        n_vars++;
      default:
        legend[Polariz] = "P_x";
        n_vars++;
    }
  }


  int PDens = -1;
  if (variables.find("PowerDensity")!= varend)
  {
    PDens = n_vars;
    legend[n_vars] = "power_density[W*cm^-3]";
    n_vars++;
  }


  legend.resize(n_vars);

  results.resize(nn * n_vars);

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");
  
  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
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

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL); 

    sc->reinit(elem);

    fe->reinit(elem);
    
    //Get the thermoelectric power------------
    double Pn =  get_electrons_thermoelectric_power(elem);
    double Pp =  get_holes_thermoelectric_power(elem);
      
    //Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_node();


    unsigned int n_dofs = dof_indices_u.size();
    // get the solution values at the centroid
    double en_x = 0.0, ep_x = 0.0;
    double en_y = 0.0, ep_y = 0.0;
    double en_z = 0.0, ep_z = 0.0;
    double dT_x = 0.0;
    double dT_y = 0.0;
    double dT_z = 0.0;
    double u  = 0.0;
    double en = 0.0;
    double ep = 0.0;
    RealGradient e_field(0);
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      en_x  += dphi[i][0](0) * solution(dof_indices_en[i]);
      en_y  += dphi[i][0](1) * solution(dof_indices_en[i]);
      en_z  += dphi[i][0](2) * solution(dof_indices_en[i]);
      
      ep_x  += dphi[i][0](0) * solution(dof_indices_ep[i]);
      ep_y  += dphi[i][0](1) * solution(dof_indices_ep[i]);
      ep_z  += dphi[i][0](2) * solution(dof_indices_ep[i]);
      
      dT_x  += dphi[i][0](0) * T_nodes[i];
      dT_y  += dphi[i][0](1) * T_nodes[i];
      dT_z  += dphi[i][0](2) * T_nodes[i];

      u  += phi[i][0] * solution(dof_indices_u[i]);
      en += phi[i][0] * solution(dof_indices_en[i]);
      ep += phi[i][0] * solution(dof_indices_ep[i]);

      e_field += dphi[i][0] * solution(dof_indices_u[i]);
    }
    e_field *= -phi0;
    en_x *= phi0;
    en_y *= phi0;
    en_z *= phi0;
    ep_x *= phi0;
    ep_y *= phi0;
    ep_z *= phi0;

    // prepare for calculating local properties
    sc->set_coordinates(elem->centroid());


    sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

    sc->set_electric_field(e_field);

    sc->calculate_densities();
    sc->calculate_mobilities();


    // we put the minus here for convenience
    double sigma_e = -Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = -Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();

    unsigned int id = n_vars * elem_number;

    if (EField != -1)
    {
      switch (dim)
      {
        case 3:
          results[id + EField + 2] = e_field(2);
        case 2:
          results[id + EField + 1] = e_field(1);
          results[id + EField + dim] = e_field.size();
        default:
          results[id + EField] = e_field(0);
      }
    }


    if (Polariz != -1)
    {
      const RealVectorValue& pol = sc->get_total_polarization();
      switch (dim)
      {
        case 3:
          results[id + Polariz + 2] = pol(2);
        case 2:
          results[id + Polariz + 1] = pol(1);
          results[id + Polariz + dim] = pol.size();
        default:
          results[id + Polariz] = pol(0);
      }
    }


    if (Jn != -1)
    {
      double jx = sigma_e * (en_x + Pn * dT_x);
      double jy = sigma_e * (en_y + Pn * dT_y);
      double jz = sigma_e * (en_z + Pn * dT_z);
      switch (dim)
      {
        case 3:
          results[id + Jn + 2] = jz;
        case 2:
          results[id + Jn + 1] = jy;
          results[id + Jn + dim] = sqrt(jx * jx + jy * jy + jz * jz);
        default:
          results[id + Jn] = jx;
      }
    }


    if (Jp != -1)
    {
      double jx = sigma_h * (ep_x + Pp * dT_x);
      double jy = sigma_h * (ep_y + Pp * dT_y);
      double jz = sigma_h * (ep_z + Pp * dT_z);
      switch (dim)
      {
        case 3:
          results[id + Jp + 2] = jz;
        case 2:
          results[id + Jp + 1] = jy;
          results[id + Jp + dim] = sqrt(jx * jx + jy * jy + jz * jz);
        default:
          results[id + Jp] = jx;
      }
    }


    if (J != -1)
    {
      double jx = sigma_e * (en_x + Pn * dT_x) + sigma_h * (ep_x + Pp * dT_x);
      double jy = sigma_e * (en_y + Pn * dT_y) + sigma_h * (ep_y + Pp * dT_y);
      double jz = sigma_e * (en_z + Pn * dT_z) + sigma_h * (ep_z + Pp * dT_z);
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


    if (PDens != -1)
    {
      double jx = sigma_e * (en_x + Pn * dT_x) + sigma_h * (ep_x + Pp * dT_x);
      double jy = sigma_e * (en_y + Pn * dT_y) + sigma_h * (ep_y + Pp * dT_y);
      double jz = sigma_e * (en_z + Pn * dT_z) + sigma_h * (ep_z + Pp * dT_z);
      double P = jx * e_field(0) + jy * e_field(1) + jz  * e_field(2);
      results[id + PDens] = P;
    }


    elem_number++;
  }

  results.resize(elem_number * n_vars);


  



}



void
DriftDiffusion::build_integrated_quantities(const set<string>& names,
    vector<double>& values)
{
  const set<string>::const_iterator varend(names.end());

  if (names.find("current") != varend)
  {
    if (get_options().current_calculation == RSTF)
      calculate_currents_rstf();
    else
      calculate_currents_surfint();

    values.resize(_boundary_currents.size());

    ContactData::iterator it(_boundary_currents.begin());
    const ContactData::iterator end(_boundary_currents.end());
    for (unsigned int id = 0; it != end; ++it, id++)
      values[id] = it->second * it->first->get_area_factor();

  }
}


void
DriftDiffusion::calculate_currents(void)
{
  if (get_options().current_calculation == RSTF)
    calculate_currents_rstf();
  else
    calculate_currents_surfint();
}


void
DriftDiffusion::build_integrated_quantities_description(
    const std::set<std::string>& names,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  const set<string>::const_iterator varend(names.end());

  if (names.find("current") != varend)
  {
    legend.resize(_boundary_currents.size());

    ContactData::iterator it(_boundary_currents.begin());
    const ContactData::iterator end(_boundary_currents.end());
    for (unsigned int id = 0; it != end; ++it, id++)
      legend[id] = it->first->get_name();

    description.resize(1);
    unsigned int dim = get_mesh().mesh_dimension();
    ostringstream s;
    s << "Contact currents. Units A";
    switch (dim)
    {
      case 1:
        s << "cm^-2";
        break;
      case 2:
        s << "cm^-1";
        break;
    }
    description[0] = s.str();
  }
}




double
DriftDiffusion::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_potential_scaling();
}



void
DriftDiffusion::set_dirichlet_bc(void)
{
  
  const Device& device = *_device;

  EquationSystems& es = get_equation_systems();

  // references for nicer code
  const Mesh& mesh = get_mesh();

  TiberNonlinearSystem& system = es.get_system<TiberNonlinearSystem>(
      get_equation_system_name());


  // the current solutions
  NumericVector<Number>& solution = system.get_solution_vector();

  BoundaryNodeList& dirichlet_nodes = _dirichlet_nodes;

  BoundaryNodeList::const_iterator node_it;
  const BoundaryNodeList::const_iterator end =
    dirichlet_nodes.end();
    

  const Scaling& scaling = get_scaling();
  const double phi0 = scaling.get_potential_scaling();


  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("potential");
  const unsigned int n_var = system.variable_number("fermi_e");
  const unsigned int p_var = system.variable_number("fermi_h");

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

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
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_n, n_var);
    dof_map.dof_indices(elem, dof_indices_p, p_var);

    unsigned int n_dofs = dof_indices_u.size();

    
    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));
    assert(sc != NULL);

    sc->reinit(elem);

    {
      // loop over all nodes and check if it is a dirichlet type node
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        node_it = dirichlet_nodes.find(elem->get_node(i));
        if (node_it != end)
        {
          Boundary* bd = node_it->second;
          ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
              bd->get_boundary_properties(get_id()));
          contact->set_material(sc);

          if (contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
          {
            double val = (contact->get_boundary_value(POTENTIAL)
                + contact->get_inner_voltage()) / phi0;
            solution.set(dof_indices_u[i], val);
          }

          if (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
          {
            double val = (contact->get_boundary_value(FERMIE)
                + contact->get_inner_voltage()) / phi0;
            solution.set(dof_indices_n[i], val);
          }

          if (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET)
          {
            double val = (contact->get_boundary_value(FERMIH)
                + contact->get_inner_voltage()) / phi0;
            solution.set(dof_indices_p[i], val);
          }
        }
      }
    }
  }
}


double DriftDiffusion::get_electron_conducibility(const Elem* elem)
{
 
  if (is_solved() && (get_environment().contains_element(elem)))   
  {  
    
    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
	 _device->get_material(elem->subdomain_id())->get_model(get_id()));


    assert(sc != NULL);


    sc->reinit(elem);


    Solution  potentials;
    get_solution(elem,elem->centroid(),potentials); 
    
    sc->set_potentials(potentials.potential,potentials.fermi_e,potentials.fermi_h);
       
    sc->calculate_densities();

    sc->calculate_mobilities();


    double sigma_e = Constants::e * sc->get_electron_density() *
        sc->get_electron_mobility();


   return sigma_e;

  }
  else
  {

    return 0.0;

  }



}

double DriftDiffusion::get_hole_conducibility(const Elem* elem)
{
  
   if (is_solved() && (get_environment().contains_element(elem)))  
  {  
    
    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
	    _device->get_material(elem->subdomain_id())->get_model(get_id()));
    
    
    assert(sc != NULL);
    
    
    sc->reinit(elem);
    

    Solution  potentials;
    get_solution(elem,elem->centroid(),potentials); 
    
    sc->set_potentials(potentials.potential,potentials.fermi_e,potentials.fermi_h);
       
    sc->calculate_densities();

    sc->calculate_mobilities();


    double sigma_h = Constants::e * sc->get_hole_density() *
        sc->get_hole_mobility();


   return sigma_h;

  }
  else
  {

    return 0.0;

  }

 


}   



double DriftDiffusion::get_electrons_thermoelectric_power(const Elem* elem)
{
     

  

  
  if (is_solved() && (get_environment().contains_element(elem))) 
  {  
    
    DriftDiffusionProperties* model =
      dynamic_cast<DriftDiffusionProperties*>(
	 _device->get_material(elem->subdomain_id())->get_model(get_id()));
    
    //initialization of the model

    double T_lattice = _device->get_temperature(elem);
  
    

    model->reinit(elem);
    
    
    //Set potentials
       
    Solution  potentials;     
    
    get_solution(elem,elem->centroid(),potentials); 
    
    model->set_potentials(potentials.potential,potentials.fermi_e,potentials.fermi_h);
    
    
    //Set Temperature
    
   
    
    model->compute_thermoelectric_powers(); 
    
    
    return (model->get_electrons_thermoelectric_power());
    
    
  }
  else
  {

    return 0.0;

    }
    
}
  




double DriftDiffusion::get_holes_thermoelectric_power(const Elem* elem)
{


  
  if (is_solved() && (get_environment().contains_element(elem))) 
  { 
    
    DriftDiffusionProperties* model =
    dynamic_cast<DriftDiffusionProperties*>(
	_device->get_material(elem->subdomain_id())->get_model(get_id()));
  
  
  //initialization of the model

  model->reinit(elem);
  
  
  //Set potentials
  
  Solution  potentials;     
  
  get_solution(elem,elem->centroid(),potentials); 
  
  model->set_potentials(potentials.potential,potentials.fermi_e,potentials.fermi_h);
  
  
  //Set Temperature
  
  
  model->compute_thermoelectric_powers(); 
  
  return (model->get_holes_thermoelectric_power());
  }
  else
  {
    
    return 0.0;
    
  }
}








  void
DriftDiffusion::assemble_system(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  switch (_this->_options.coupling)
  {
    case (POISSON | ECURRENT):
      _this->do_assembly<POISSON | ECURRENT>(x, residual, jacobian);
      break;
    case (POISSON | HCURRENT):
      _this->do_assembly<POISSON | HCURRENT>(x, residual, jacobian);
      break;
    case (CURRENTS):
      _this->do_assembly<CURRENTS>(x, residual, jacobian);
      break;
    case (POISSON):
      _this->do_assembly<POISSON>(x, residual, jacobian);
      break;
    case (ECURRENT):
      _this->do_assembly<ECURRENT>(x, residual, jacobian);
      break;
    case (HCURRENT):
      _this->do_assembly<HCURRENT>(x, residual, jacobian);
      break;
    default:
      _this->do_assembly<FULLYCOUPLED>(x, residual, jacobian);
  }

}






template <int coupling>
void
DriftDiffusion::do_assembly(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  PerfLog perf_log("Matrix assembly", false);
  perf_log.start_event("assembly");
  
  // references for nicer code
  const Mesh& mesh = get_mesh();
  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = static_cast<TiberNonlinearSystem&>(
      eq_sys.get_system(get_equation_system_name()));

  const unsigned int dim = mesh.mesh_dimension();
  
  const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  const Options& params = get_options();
  Options& options = get_options();

  BoundaryNodeList& dirichlet_nodes = _dirichlet_nodes;

  const NumericVector<Number>& dx = *(system.solution);

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
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  // x 1e4 because we calculate in cm
  const double P0 = (Constants::e * x0 * C0) * 1e4;
  // density scaling for electrons
  double C0_e = C0;
  // density scaling for holes
  double C0_h = C0;
 
  if (do_local_scaling_)
    C0_e = C0_h = 1.0;


  // scaling for recombination rates
  double R0_e = C0_e / scaling.get_time_scaling();
  double R0_h = C0_h / scaling.get_time_scaling();


  const DofMap& dof_map = system.get_dof_map();
  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_e");
  const unsigned int ep_var = system.variable_number("fermi_h");
  
  FEType fe_type = system.variable_type(u_var);

  libMeshEnums::Order integration_order = params.integration_order;

  // the finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, integration_order);
  fe->attach_quadrature_rule(&qrule);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  
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



  // references to boundary-specific data that will be used to
  // assemble the system.
  // Data will be given for each quadrature point.
  // 
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  //
  const vector<vector<RealGradient> >&  dphi_face = fe_face->get_dphi();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point_face = fe_face->get_xyz();
  //
  const vector<Point>& face_normals = fe_face->get_normals();
  //
  // Jacobian * quadrature weight at each integration point.   
  const vector<Real>& JxW_face = fe_face->get_JxW();


  // the system matrix (will hold also element jacobian contribution)
  DenseMatrix<Number> Ke;
  // the system rhs (will hold also element rhs contribution)
  DenseVector<Number> Fe;
  // the local solution
  DenseVector<Number> X;
  // the local old step
  DenseVector<Number> dX;

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

  DenseSubVector<Number>
    dXu(dX),
    dXn(dX),
    dXp(dX);


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

    ID subdomain = elem->subdomain_id();

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
    dX.resize(n_dofs_tot);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);
    dof_map.extract_local_vector(dx, dof_indices, dX);

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
    //
    dXu.reposition(0, n_dofs);
    dXn.reposition(n_dofs, n_dofs);
    dXp.reposition(2 * n_dofs, n_dofs);



    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);

    sc->reinit(elem);
    

    //Get the thermoelectric power
    double eTEpower =  get_electrons_thermoelectric_power(elem) / phi0;
    double hTEpower =  get_holes_thermoelectric_power(elem) / phi0;

           
     //Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_node();
   

    vector<vector<double> > local_scaling(elem->n_nodes(), vector<double>(3, 1));
    if (do_local_scaling_)
    {
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        local_scaling[n] = local_scaling_[elem->get_node(n)];
        //local_scaling[n][2] = 1.0;
      }
    }




    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      RealGradient e_field(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * Xu(i);
        en += phi[i][qp] * Xn(i);
        ep += phi[i][qp] * Xp(i);
        e_field += dphi[i][qp] * Xu(i);
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

      sc->set_electric_field(phi0 * e_field);

      sc->calculate_densities();

      long double n = sc->get_electron_density();
      long double p = sc->get_hole_density();
      //double Nd = sc->get_ionized_donor_density();
      //double Na = sc->get_ionized_acceptor_density();
  
      // calculate all local properties
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

    
      double epsilon = sc->get_relative_permittivity();
      long double l2_eps = l2 * epsilon;

      long double Rn = sc->get_net_electron_recombination_rate();
      //Rn = (fabs(Rn) < 1.0e-19) ? 0.0 : Rn;
      long double Rp = sc->get_net_hole_recombination_rate();
      //Rp = (fabs(Rp) < 1.0e-19) ? 0.0 : Rp;
      

      //double ni = sc->get_intrinsic_density();
      long double mue = sc->get_electron_mobility();
      long double muh = sc->get_hole_mobility();


      // the jacobian x weight x scaling
      long double J = JxW[qp];


      // NOTE: sigma_e = mu_e * n is the electron conductivity
      long double sigma_e = mue * n / (mu0 * C0_e);
      long double sigma_h = muh * p / (mu0 * C0_h);
      long double sigma_e_x_Pe_x_J = J * sigma_e * eTEpower;
      long double sigma_h_x_Ph_x_J = J * sigma_h * hTEpower;

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

      // 
      // First we will build the system matrix Ke_ij
      //
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          long double laplace =
            J * (dphi[i][qp] * dphi[j][qp]);
          
          if (coupling & POISSON)
            Kuu(i,j) += l2_eps * laplace / local_scaling[i][2];
          
          if (coupling & ECURRENT)
            Knn(i,j) += sigma_e * laplace / local_scaling[i][0];
          
          if (coupling & HCURRENT)
            Kpp(i,j) += sigma_h * laplace / local_scaling[i][1];
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
        long double dn_dphi = sc->get_electron_density_derivative();
        long double dp_dphi = sc->get_hole_density_derivative();
        long double dNd_dphi = sc->get_ionized_donor_density_derivative();
        long double dNa_dphi = sc->get_ionized_acceptor_density_derivative();

        long double drho[3];
        drho[1] = (dn_dphi - dNd_dphi) * phi0 / C0;
        drho[2] = -(dp_dphi - dNa_dphi) * phi0 / C0;
        drho[0] = -(drho[1] + drho[2]);
        if (sc->is_dielectric())
          drho[2] = drho[1] = drho[0] = 0.0;

        long double dRn_dn =
          sc->get_net_electron_recombination_rate_derivatives()[0];
        long double dRn_dp =
          sc->get_net_electron_recombination_rate_derivatives()[1];
        long double dRp_dn = sc->get_net_hole_recombination_rate_derivatives()[0];
        long double dRp_dp = sc->get_net_hole_recombination_rate_derivatives()[1];

        long double dRn[3];
        long double dRp[3];
        dRn[1] = -dRn_dn * dn_dphi * phi0 / R0_e;
        dRn[2] = -dRn_dp * dp_dphi * phi0 / R0_e;
        dRn[0] = -(dRn[1] + dRn[2]);
        dRp[1] = -dRp_dn * dn_dphi * phi0 / R0_h;
        dRp[2] = -dRp_dp * dp_dphi * phi0 / R0_h;
        dRp[0] = -(dRp[1] + dRp[2]);

        if (Rn == 0.0)
          dRn[0] = dRn[1] = dRn[2] = 0.0;
        if (Rp == 0.0)
          dRp[0] = dRp[1] = dRp[2] = 0.0;

        // d(sigma_n)/du * element-jacobian
        // sigma_n = mu_n * n means the conductivity of electrons
        long double dsigma_e = J * phi0 / (mu0 * C0_e) * mue * dn_dphi;
        long double dsigma_h = J * phi0 / (mu0 * C0_h) * muh * dp_dphi;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)

            if (_options.exact_newton)
            {
              long double dsigma_e_x_phi = dsigma_e * phi[j][qp];
              long double dsigma_h_x_phi = dsigma_h * phi[j][qp];

              for (unsigned int k = 0; k < n_dofs; k++)
              {
                long double laplace = (dphi[i][qp] * dphi[k][qp]);

                if (coupling & ECURRENT)
                {
                  long double elem_contrib =
                    dsigma_e_x_phi * laplace * Xn(k) / local_scaling[i][0];

                  if (coupling & POISSON)
                    Knu(i,j) += elem_contrib;

                  Knn(i,j) -= elem_contrib;
                }

                if (coupling & HCURRENT)
                {
                  long double elem_contrib =
                    dsigma_h_x_phi * laplace * Xp(k) / local_scaling[i][1];

                  if (coupling & POISSON)
                    Kpu(i,j) += elem_contrib;

                  Kpp(i,j) -= elem_contrib;
                }
              }
            }


            // contribution of the Seebeck effect ->residual_derivative
            double dsigma_e_x_phi_x_Pe = dsigma_e * phi[j][qp] * eTEpower;
            double dsigma_h_x_phi_x_Ph = dsigma_h * phi[j][qp] * hTEpower;
	    
            for (unsigned int k = 0; k < n_dofs; k++)
            {
              double laplace = dphi[i][qp] * dphi[k][qp];
	      
              if (coupling & ECURRENT)
              {
                double elem_contrib =
                  dsigma_e_x_phi_x_Pe * laplace * T_nodes[k] / local_scaling[i][0];
		
		Knn(i,j) -= elem_contrib;
		
                if (coupling & POISSON)
		  Knu(i,j) += elem_contrib;
		
              }
	      
              if (coupling & HCURRENT)
              {
                double elem_contrib =
                  dsigma_h_x_phi_x_Ph * laplace * T_nodes[k] / local_scaling[i][1];
		
		Kpp(i,j) -= elem_contrib;
		
		if (coupling & POISSON)
		  Kpu(i,j) += elem_contrib;
					     
		}
            }
	    




            // The dFe_i/dX_j part
            long double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (_options.exact_newton || (j == i))
            {
              if (coupling & POISSON)
              {
                Kuu(i,j) -= drho[0] * phi_i_x_phi_j / local_scaling[i][2];

                if (coupling & ECURRENT)
                  Kun(i,j) -= drho[1] * phi_i_x_phi_j / local_scaling[i][2];

                if (coupling & HCURRENT)
                  Kup(i,j) -= drho[2] * phi_i_x_phi_j / local_scaling[i][2];
              }            

              if (coupling & ECURRENT)
              {
                // could destroy M-Matrix property
                if (coupling & POISSON)
                  Knu(i,j) -= dRn[0] * phi_i_x_phi_j / local_scaling[i][0];

                Knn(i,j) -= dRn[1] * phi_i_x_phi_j / local_scaling[i][0];

                if (coupling & HCURRENT)
                  Knp(i,j) -= dRn[2] * phi_i_x_phi_j / local_scaling[i][0];
              }

              if (coupling & HCURRENT)
              {
                // could destroy M-Matrix property
                if (coupling & POISSON)
                  Kpu(i,j) += dRp[0] * phi_i_x_phi_j / local_scaling[i][1];

                if (coupling & ECURRENT)
                  Kpn(i,j) += dRp[1] * phi_i_x_phi_j / local_scaling[i][1];

                Kpp(i,j) += dRp[2] * phi_i_x_phi_j / local_scaling[i][1];
              }
            }

          }
        }
/*
        if (residual != NULL)
        {
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
            {

              double dsigma_e_x_phi = dsigma_e * phi[j][qp] * 0.1;
              double dsigma_h_x_phi = dsigma_h * phi[j][qp] * 0.1;


              for (unsigned int k = 0; k < n_dofs; k++)
              {
                double laplace = (dphi[i][qp] * dphi[k][qp]);

                if (coupling & ECURRENT)
                {
                  double elem_contrib =
                    dsigma_e_x_phi * laplace * Xn(k);

                  if (coupling & POISSON)
                    Fn(i) += elem_contrib * dXu(j);

                  Fn(i) -= elem_contrib * dXn(j);
                }

                if (coupling & HCURRENT)
                {
                  double elem_contrib =
                    dsigma_h_x_phi * laplace * Xp(k);

                  if (coupling & POISSON)
                    Fp(i) += elem_contrib * dXu(j);

                  Fp(i) -= elem_contrib * dXp(j);
                }
              }


              if (i != j)
              {
                double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp] * 0.1;

                if (coupling & POISSON)
                {
                  Fu(i) -= drho[0] * phi_i_x_phi_j * dXu(j);

                  if (coupling & ECURRENT)
                    Fu(i) -= drho[1] * phi_i_x_phi_j * dXn(j);

                  if (coupling & HCURRENT)
                    Fu(i) -= drho[2] * phi_i_x_phi_j * dXp(j);
                }            

                if (coupling & ECURRENT)
                {
                  // (1) would destroy M-Matrix property
                  // (2) is 0 for Boltzmann statistics
                  //if (coupling & POISSON)
                  //  Knu(i,j) -= dRn[0] * phi_i_x_phi_j;

                  Fn(i) -= dRn[1] * phi_i_x_phi_j / local_scaling[i][0] * dXn(j);

                  if (coupling & HCURRENT)
                    Fn(i) -= dRn[2] * phi_i_x_phi_j / local_scaling[i][0] * dXp(j);
                }

                if (coupling & HCURRENT)
                {
                  // (1) would destroy M-Matrix property
                  // (2) is 0 for Boltzmann statistics
                  //if (coupling & POISSON)
                  //  Kpu(i,j) += dRp[0] * phi_i_x_phi_j;

                  if (coupling & ECURRENT)
                    Fp(i) += dRp[1] * phi_i_x_phi_j / local_scaling[i][1] * dXn(j);

                  Fp(i) += dRp[2] * phi_i_x_phi_j / local_scaling[i][1] * dXp(j);
                }
              }
            }
          }
        } // if (residual != NULL)
*/

      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        long double J_x_rho = J * sc->get_charge_density() / C0;
        if (sc->is_dielectric())
          J_x_rho = 0.0;

        long double J_x_P0 = J / P0;

        // net recombination rate
        long double J_x_Rn = J * Rn / R0_e;
        long double J_x_Rp = J * Rp / R0_h;

        RealVectorValue P(sc->get_total_polarization());
        P *= J_x_P0;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          long double net_recomb_e = J_x_Rn * phi[i][qp] / local_scaling[i][0];
          long double net_recomb_h = J_x_Rp * phi[i][qp] / local_scaling[i][1];
          
          if (coupling & POISSON)
            Fu(i) -= (J_x_rho * phi[i][qp] + (P * dphi[i][qp]))
              / local_scaling[i][2];
          else
            Fu(i) -= Xu(i);
          
          if (coupling & ECURRENT)
            Fn(i) -= net_recomb_e;
          else
            Fn(i) -= Xn(i);

          if (coupling & HCURRENT)
            Fp(i) += net_recomb_h;
          else
            Fp(i) -= Xp(i);
        }
	


        // include Seebeck contribution -> Residual
	for (unsigned int i = 0; i < n_dofs; i++)
	{
	  for (unsigned int k = 0; k < n_dofs; k++)
	  {
	    Real laplace = dphi[i][qp] * dphi[k][qp];
	    
	    if (coupling & ECURRENT)
	      Fn(i) += sigma_e_x_Pe_x_J * laplace *
		T_nodes[k] / local_scaling[i][0];
	    
	    if (coupling & HCURRENT)
	      Fp(i) += sigma_h_x_Ph_x_J * laplace *
		T_nodes[k] / local_scaling[i][1];
	  }
	}
	
	
	
	
      }
    } // end loop over quadrature points
    
    
    
    // this is for surface resistance
    map<unsigned int, double> nodal_flux_n;
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
      nodal_flux_n[elem->node(i)] = 0.0;
    map<unsigned int, double> nodal_flux_p(nodal_flux_n);

    set<unsigned int> nodes_on_boundary_sides;

    
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
      ElementSide side(top_parent, s);
      
      const Elem* neighbour = elem->neighbor(s);

      // is this a boundary?
      if (environment.is_boundary(side))
      {
        // we need to know if it is an outer boundary
        bool true_boundary = environment.is_outer_boundary(side);

        Boundary* boundary = environment.get_boundary(side);

        ElectricalContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<ElectricalContact*>(
              boundary->get_boundary_properties(get_id()));


        fe_face->reinit(elem, s);

        // calculate the fluxes on the nodes
        if ((contact != NULL) && (dim > 1))
        {
          AutoPtr<Elem> side(elem->build_side(s));
          
          vector<Point> p(side->n_nodes());
          for (unsigned int i = 0; i < side->n_nodes(); i++)
          {
            nodes_on_boundary_sides.insert(side->node(i));
            p[i] = side->point(i);
          }
          fe->reinit(elem, &p);

          for (unsigned int i = 0; i < side->n_nodes(); i++)
          {
            RealGradient e_field(0);
            RealGradient grad_en(0);
            RealGradient grad_ep(0);
            RealGradient grad_T(0);
            for (unsigned int j = 0; j < n_dofs; j++)
            {
              e_field += dphi[j][i] * Xu(j);
              grad_en += dphi[j][i] * Xn(j);
              grad_ep += dphi[j][i] * Xp(j);
            }

            sc->set_potentials(phi0 * Xu(i), phi0 * Xn(i), phi0 * Xp(i));
            sc->set_coordinates(side->point(i));
            sc->set_electric_field(-phi0 / x0 * e_field);
            sc->calculate_densities();
            sc->calculate_mobilities();

            // we put phi0 here for convenience
            double sigma_e = phi0 * sc->get_electron_mobility() *
              sc->get_electron_density();
            double sigma_h = phi0 * sc->get_hole_mobility() *
              sc->get_hole_density();

            double Pn =  get_electrons_thermoelectric_power(elem);
            double Pp =  get_holes_thermoelectric_power(elem);

            if (coupling & ECURRENT)
              nodal_flux_n[side->node(i)] = 
                (sigma_e * grad_en + Pn * grad_T) * face_normals[0] / x0;
            if (coupling & HCURRENT)
              nodal_flux_p[side->node(i)] = 
                -(sigma_h * grad_ep + Pp * grad_T) * face_normals[0] / x0;
          }
        }


        // for von Neumann or mixed type boundary conditions
        vector<double> coeff(3, 0.0);
        vector<double> value(3, 0.0);

        // the derivatives
        vector<vector<double> > dcoeff(3);
        vector<vector<double> > dvalue(3);

        //
        // NOTE: we have to integrate over the boundary also if there are
        //       no contacts because there could be polarization.
        //

        if (dim > 1)
        {

          int phi_size = phi_face.size();

          // now integrate to include von Neumann and mixed type BCs
          // and polarization
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
   
            double epsilon = sc->get_relative_permittivity();
            double l2_eps = l2 * epsilon;

            // get the boundary condition coefficients
            if (contact != NULL)
            {

              // get the solution values at the quadrature point
              Real u  = 0.0;
              Real en = 0.0;
              Real ep = 0.0;
              RealGradient e_field(0);
              RealGradient grad_en(0);
              RealGradient grad_ep(0);
              RealGradient grad_T(0);
              for (unsigned int i = 0; i < n_dofs; i++)
              {
                u  += phi_face[i][qp] * Xu(i);
                en += phi_face[i][qp] * Xn(i);
                ep += phi_face[i][qp] * Xp(i);
                e_field += dphi_face[i][qp] * Xu(i);
                grad_en += dphi_face[i][qp] * Xn(i);
                grad_ep += dphi_face[i][qp] * Xp(i);
              }

              sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
              sc->set_coordinates(q_point_face[qp]);
              sc->set_electric_field(phi0 / x0 * e_field);
              sc->calculate_densities();
              sc->calculate_mobilities();

              // we put the phi0 here for convenience
              double sigma_e = phi0 * sc->get_electron_mobility() *
                sc->get_electron_density();
              double sigma_h = phi0 * sc->get_hole_mobility() *
                sc->get_hole_density();

              double Pn =  get_electrons_thermoelectric_power(elem);
              double Pp =  get_holes_thermoelectric_power(elem);

              double jn = 0.0;
              double jp = 0.0;
              if (coupling & ECURRENT)
                jn = (sigma_e * grad_en + Pn * grad_T) * face_normals[qp] / x0;
              if (coupling & HCURRENT)
                jp = -(sigma_h * grad_ep + Pp * grad_T) * face_normals[qp] / x0;


              contact->set_material(sc);
              contact->set_normal_fluxes(jn, jp);
              double a, c;

              if (coupling & POISSON)
              {
                contact->get_normal_derivative(POTENTIAL, a, c);
                contact->get_derivatives_of_normal_derivative(POTENTIAL,
                    dcoeff[0], dvalue[0]);
                //coeff[0] = a * x0;
                coeff[0] = 0.0;
                //value[0] = c * x0 / phi0;
                value[0] = c / (x0 * C0);
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
            double J = JxW_face[qp];

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
/*
              //double val_uu = J * (dcoeff[0][0] * u - dvalue[0][0]);
              double val_uu = J * (- dvalue[0][0] / x0 / C0 * phi0);

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                for (unsigned int j = 0; j < n_dofs; j++)
                {

                  Real phi_i_x_phi_j =
                    phi_face[i][qp] * phi_face[j][qp];

                  if (coupling & POISSON)
                    Kuu(i,j) += coeff[0] * phi_i_x_phi_j;

                  if (coupling & ECURRENT)
                    Knn(i,j) += coeff[1] * phi_i_x_phi_j;

                  if (coupling & HCURRENT)
                    Kpp(i,j) += coeff[2] * phi_i_x_phi_j;
                }
              }
*/
            }

            // contribution to -Fe_i
            if (residual != NULL)
            {
              double Pn = 0.0;
              if (true_boundary && (contact == NULL))
              {
                // If we are on an outer boundary, we have to include
                // the polarization
                //
                // NOTE:
                // we only include the polarization when no boundary
                // is defined
                RealVectorValue P(sc->get_total_polarization());
                Pn = (P * face_normals[qp]) / P0;
              }
              //double value_u = J * (l2_eps * value[0] - Pn);
              double value_u = J * (value[0] - Pn);
              double value_n = J * value[1] / (mu0 * C0_e);
              double value_p = J * value[2] / (mu0 * C0_h);

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                if (coupling & POISSON)
                  Fu(i) -= value_u * phi_face[i][qp] / local_scaling[i][2];

                if (coupling & ECURRENT)
                  Fn(i) -= value_n * phi_face[i][qp] / local_scaling[i][0];

                if (coupling & HCURRENT)
                  Fp(i) -= value_p * phi_face[i][qp] / local_scaling[i][1];
              }
            } 
          }
        }
        else // i.e. dim == 1
        {
          nodes_on_boundary_sides.insert(elem->node(s));

          // s is the node of the element lying on the boundary
          Real u  = Xu(s);
          Real en = Xn(s);
          Real ep = Xp(s);

          // calculate densities etc.
          sc->set_coordinates(elem->point(s));
          sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
          sc->set_coordinates(elem->point(s));
          
          RealGradient e_field(0.0);
          double grad_en = 0.0;
          double grad_ep = 0.0;
          double grad_T = 0.0;
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            grad_en += dphi_face[n][0](0)* Xn(n);
            grad_ep += dphi_face[n][0](0) * Xp(n);
            e_field(0) += dphi_face[n][0](0) * Xu(n);

            grad_T += dphi_face[n][0](0) * T_nodes[n];
          }
          
          sc->set_electric_field(phi0 / x0 * e_field);
          sc->calculate_densities();
          sc->calculate_mobilities();

          // we put the phi0 here for convenience
          double sigma_e = phi0 * sc->get_electron_mobility() *
            sc->get_electron_density();
          double sigma_h = phi0 * sc->get_hole_mobility() *
            sc->get_hole_density();

          double Pn =  get_electrons_thermoelectric_power(elem);
          double Pp =  get_holes_thermoelectric_power(elem);

          double x_c = elem->centroid()(0);
          double x_s = elem->point(s)(0);
          double sign = (x_s > x_c) ? 1 : -1;

          double jn = sign * (sigma_e * grad_en + Pn * grad_T) / x0;
          double jp = -sign * (sigma_h * grad_ep + Pp * grad_T) / x0;
          
          double epsilon = sc->get_relative_permittivity();
          double l2_eps = l2 * epsilon;

          if (coupling & ECURRENT)
            nodal_flux_n[elem->node(s)] = jn;
          if (coupling & HCURRENT)
            nodal_flux_p[elem->node(s)] = jp;

          // get the boundary condition coefficients
          if (contact != NULL)
          {
            contact->set_material(sc);
            contact->set_normal_fluxes(jn, jp);

            double a, c;

            if (coupling & POISSON)
            {
              contact->get_normal_derivative(POTENTIAL, a, c);
              contact->get_derivatives_of_normal_derivative(POTENTIAL,
                  dcoeff[0], dvalue[0]);
              //coeff[0] = a * x0;
              coeff[0] = 0.0;
              //value[0] = c * x0 / phi0;
              value[0] = c / (x0 * C0);
            }
            if (coupling & ECURRENT)
            {
              contact->get_normal_derivative(FERMIE, a, c);
              contact->get_derivatives_of_normal_derivative(FERMIE,
                  dcoeff[1], dvalue[1]);
              coeff[1] = a * x0;
              value[1] = c * x0 / phi0;
            }
            if (coupling & HCURRENT)
            {
              contact->get_normal_derivative(FERMIH, a, c);
              contact->get_derivatives_of_normal_derivative(FERMIH,
                  dcoeff[2], dvalue[2]);
              coeff[2] = a * x0;
              value[2] = c * x0 / phi0;
            }
          }


          // first the contributions to Ke_ij
          if (coupling & POISSON)
            Kuu(s,s) += l2_eps * coeff[0] / local_scaling[s][2];

          if (coupling & ECURRENT)
            Knn(s,s) += coeff[1] / local_scaling[s][0];

          if (coupling & HCURRENT)
            Kpp(s,s) += coeff[2] / local_scaling[s][1];


          // contribution to the jacobian
          if (jacobian != NULL)
          {
            double val_uu = J * (- dvalue[0][0] / x0 / C0 * phi0);
            double val_nn = J * (- dvalue[1][1] / mu0 / C0 * phi0);
            double val_pp = J * (- dvalue[2][2] / mu0 / C0 * phi0);

            if (coupling & POISSON)
              Kuu(s,s) += val_uu / local_scaling[s][2];

            if (coupling & ECURRENT)
              Knn(s,s) += val_nn / local_scaling[s][0];

            if (coupling & HCURRENT)
              Kpp(s,s) += val_pp / local_scaling[s][1];
          }

          
          // contribution to -Fe_i
          if (residual != NULL)
          {
            double Pn =  0.0;
            if (true_boundary && (contact == NULL))
            {
              // If we are on an outer boundary, we have to include
              // the polarization
              //
              // NOTE:
              // we only include the polarization when no boundary
              // is defined
              Pn = sc->get_total_polarization()(0) / P0;

              // what is the outer normal in this point??
              // Idea: if x(s) > x(centroid), normal is +1
              //       else it is -1
              double x_c = elem->centroid()(0);
              double x_s = elem->point(s)(0);
              Pn = (x_s > x_c) ? Pn : -Pn;
            }
            //double value_u = l2_eps * value[0] - Pn;
            double value_u = value[0] - Pn;
            double value_n = value[1] / (mu0 * C0_e);
            double value_p = value[2] / (mu0 * C0_h);


            if (coupling & POISSON)
              Fu(s) -= value_u / local_scaling[s][2];

            if (coupling & ECURRENT)
              Fn(s) -= value_n / local_scaling[s][0];

            if (coupling & HCURRENT)
              Fp(s) -= value_p / local_scaling[s][1];
          }
        }
      }
    } // end loop over element sides



    // check if it is a dielectric
    if (sc->is_dielectric())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          Kun(i, j) = Kup(i, j) = 0.0;
          Knu(i, j) = Knn(i, j) = Knp(i, j) = 0.0;
          Kpu(i, j) = Kpn(i, j) = Kpp(i, j) = 0.0;
        }

        if (is_dielectric_boundary_node(elem->get_node(i)))
          Knn(i, i) = Kpp(i, i) = 0.0;
        else
          Knn(i, i) = Kpp(i, i) = 1.0;

        // we simply set the electrochemical potential to zero
        Fn(i) = 0.0;
        Fp(i) = 0.0;
      }
    }


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
          Boundary* bd = node_it->second;
          ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
              bd->get_boundary_properties(get_id()));
          contact->set_material(sc);
          contact->set_normal_fluxes(
              nodal_flux_n[elem->node(i)], nodal_flux_p[elem->node(i)]);

          // we only impose Dirichlet type BCs if the node has an associated
          // boundary side
          if (nodes_on_boundary_sides.find(elem->node(i)) !=
              nodes_on_boundary_sides.end())
          {

            if (coupling & POISSON)
            {
              if (contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
              {
                double val = (contact->get_boundary_value(POTENTIAL)
                    + contact->get_inner_voltage()) / phi0;
                Ke.condense(i, i, -val, Fe);
              }
            }

            if (coupling & ECURRENT)
            {
              if (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
              {
                double val = (contact->get_boundary_value(FERMIE)
                    + contact->get_inner_voltage()) / phi0;
                Ke.condense(i + n_dofs, i + n_dofs, -val, Fe);
              }
            }

            if (coupling & HCURRENT)
            {
              if (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET)
              {
                double val = (contact->get_boundary_value(FERMIH)
                    + contact->get_inner_voltage()) / phi0;
                Ke.condense(i + 2 * n_dofs, i + 2 * n_dofs, -val, Fe);
              }
            }
          }
          else
          {
            // in this case we do not change the matrix and rhs

            if (coupling & POISSON)
            {
              if (contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
              {
                for (int j = 0; j < n_dofs_tot; j++)
                  Ke(i, j) = 0.0;
                Fe(i) = 0;
              }
            }

            if (coupling & ECURRENT)
            {
              if (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
              {
                for (int j = 0; j < n_dofs_tot; j++)
                  Ke(i + n_dofs, j) = 0.0;
                Fe(i + n_dofs) = 0;
              }
            }

            if (coupling & HCURRENT)
            {
              if (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET)
              {
                for (int j = 0; j < n_dofs_tot; j++)
                  Ke(i + 2 * n_dofs, j) = 0.0;
                Fe(i + 2 * n_dofs) = 0;
              }
            }

          }
        }
      }
    }
    else
    {
      // TODO this needs to be checked!!!
      // TODO Is now probably broken

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
            Boundary* bd = node_it->second;
            ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
                bd->get_boundary_properties(get_id()));
            contact->set_material(sc);
            contact->set_normal_fluxes(0.0, 0.0);

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
                        + contact->get_inner_voltage()) / phi0;
                    Ke.condense(id, id, -val, Fe);
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
                        + contact->get_inner_voltage()) / phi0;
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
                        + contact->get_inner_voltage()) / phi0;
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

  if (jacobian != NULL)
    jacobian->close();
  else
    residual->close();

/*
  if (coupling & ECURRENT)
  //if (1)
  {
    static int cnt = 0;
    static int cnt2 = 0;
    if (jacobian != NULL)
    {
      cerr << "write jacobian." << endl;
      ostringstream s;
      s << "jac_" << cnt << ".m";
      jacobian->print_matlab(s.str());
      cnt++;
    }
    else
    {
      cerr << "write residual." << endl;
      ostringstream s;
      s << "res_" << cnt2 << ".m";
      residual->print_matlab(s.str());
      ostringstream so;
      so << "sol_" << cnt2 << ".m";
      x.print_matlab(so.str());
      cnt2++;
    }
  }
*/  
  perf_log.stop_event("assembly");
} 








//
// explicit instantiations of template methods
//

template void
DriftDiffusion::get_solution<double>(const Elem* elem, const Point& p,
    double& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::Solution>(const Elem* elem,
    const Point& p, DriftDiffusion::Solution& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::Currents>(const Elem* elem,
    const Point& p, DriftDiffusion::Currents& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::EField>(const Elem* elem,
    const Point& p, DriftDiffusion::EField& solution);


template void
DriftDiffusion::get_solution<double>(const Elem* elem, const vector<Point>& p,
    vector<double>& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::Solution>(const Elem* elem,
    const vector<Point>& p, vector<DriftDiffusion::Solution>& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::Currents>(const Elem* elem,
    const vector<Point>& p, vector<DriftDiffusion::Currents>& solution);

template void
DriftDiffusion::get_solution<DriftDiffusion::EField>(const Elem* elem,
    const vector<Point>& p, vector<DriftDiffusion::EField>& solution);




