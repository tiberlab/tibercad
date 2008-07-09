// $Id$

// module includes
#include "DSSC.h"
#include "SimulationEnvironment.h"
#include "Scaling.h"
#include "Material.h"
#include "Boundary.h"
#include "Constants.h"
#include "DSSCModel.h"
#include "DSSCContact.h"
#include "TiberNonlinearSystem.h"
#include "SolveFailedException.h"


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


//
// Module interface
//

TIBER_MODULE(DSSC,dssc)




DSSC*
DSSC::_this;





DSSC::DSSC(void)
  : _rebuild_eq_system(true),
    _poisson_only(false),
    _scaling_type(Scaling::UNITS)
{
}




DSSC::~DSSC(void)
{
  cleanup_solver();
}
 



PhysicalModel*
DSSC::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  string modelname;
  
  modelname = options.get_option("model", "default");
  
  DSSCModel* model =
    DSSCModel::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "DSSC: No such physical model: " + modelname);

  return model;
}





BoundaryProperties*
DSSC::create_boundary_model(const ModelOptions& options) const
throw (ModelErrorException)
{
  const string& modelname = options.get_option("type", "");

  BoundaryProperties* model = NULL;

  //if (modelname != "")
  //{
    model =  DSSCContact::create(modelname, options);

    if (model == NULL)
      throw ModelErrorException(
          "DSSC: No such boundary model: " + modelname);
  //}

  return model;
}




void
DSSC::compute_scaling(Scaling::ScalingType type)
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
  double mu0 = 1;
  double C0 = 1;
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
    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          _device->get_material(elem->subdomain_id())->get_model(get_id()));

    sc->set_coordinates(elem->centroid());
    sc->set_potentials(0.0);
    sc->set_electric_field(RealGradient(0));
    sc->set_grad_fermi_n(RealGradient(0));
    sc->set_grad_fermi_I(RealGradient(0));
    sc->set_grad_fermi_I3(RealGradient(0));
    sc->set_grad_fermi_C(RealGradient(0));
    sc->reinit(elem);

    sc->calculate_densities();


    
    /*
    double mu = sc->get_hole_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;
    mu = sc->get_electron_mobility();
    mu0 = (mu0 > mu) ? mu0 : mu;
    */

    double C = sc->get_equilibrium_concentrations().n;
    C0 = (C0 > C) ? C0 : C;

    //double ni = sc->get_intrinsic_density();
    //ni0 = (ni0 > ni) ? ni0 : ni;

    double eps = sc->get_relative_permittivity();
    eps0 = (eps0 > eps) ? eps0 : eps;
  }

  switch (type)
  {
    default: // UNITS

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
  get_scaling().set_length_scaling(x0 * mesh_units);
  get_scaling().set_mobility_scaling(mu0);
  get_scaling().set_density_scaling(C0);
}






void
DSSC::reset_solver(void)
{
  if (!_rebuild_eq_system)
  {
    //get_equation_systems().delete_system(get_equation_system_name());
    _rebuild_eq_system = true;
  }
}




void
DSSC::cleanup_solver(void)
{

  // erase boundary current data structure
  _boundary_currents.erase(_boundary_currents.begin(),
      _boundary_currents.end());

  // erase dirichlet nodes data structure
  _dirichlet_nodes.erase(_dirichlet_nodes.begin(),
      _dirichlet_nodes.end());

  reset_solver();
}




void
DSSC::do_solve(void)
{
  
  cout << endl;
  cout << "<<-------------------------------------------------------------------"
    << endl;
  cout << "DSSC (name: " << get_name() << ")" << endl;

  // rebuild the system if needed
  //rebuild_equation_system();

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;
  
  parse_options();


  if (!equilibrium_done())
  {
    solve_equilibrium();
  }



  if (do_local_scaling_)
    build_local_scaling();


  
  try
  {
    do_newton();
  }
  catch (SolverException& e)
  {
    string msg = "solve failed (" +
      string(e.what()) + ")";
    throw SolveFailedException(msg);
  }

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
    /*
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
    */
  }
  cout << endl;
  cout << "------------------------------------------------------------------->>"
    << endl;

}




void
DSSC::do_equilibrium(void)
{

  // rebuild the system if needed
  //rebuild_equation_system();

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;
  
  //parse_options();


  // first we have to compute the scaling
  compute_scaling(_scaling_type);

/*
  TiberNonlinearSystem& system =
    get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());
  
  ModelOptions& solveropts = get_solver_options();
  int max_it = solveropts.get_option("nonlin_max_it", -1);
  solveropts.set_option("nonlin_max_it", 150);


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
    cout << "Solving equilibrium" << endl;

    do_newton();
    
    cout << "Equilibrium done" << endl;
  }
  catch (runtime_error& e)
  {
    cerr << "Equilibrium did not converge: " << e.what() << endl;
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

  if (max_it != -1)
    solveropts.set_option("nonlin_max_it", max_it);
  else
    solveropts.delete_option("nonlin_max_it");
*/
}



void
DSSC::guess_equilibrium(void)
{
/*
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
    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          _device->get_material(elem->subdomain_id())->get_model(get_id()));

    dof_map_u.dof_indices(elem, dof_indices_u, u_var);
    for (int i = 0; i < elem->n_nodes(); i++)
    {
      solution_u.add(dof_indices_u[i],
          sc->get_equilibrium_fermi_level()
          / (phi0 * static_cast<Real>(node_conn[elem->node(i)])));
    }
  }
*/
}



void
DSSC::do_print_info(void)
{
  parse_const_options();
  parse_options();

  string space("  ");

  if (do_local_scaling_)
    cout << space << "using local scaling" << endl;
}





void
DSSC::parse_const_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();

  do_local_scaling_ = opts.get_option("local_scaling", false);

  string scaling = opts.get_option("scaling", "");
  if (scaling == "none")
    _scaling_type = Scaling::NONE;
  else
    _scaling_type = Scaling::UNITS;


}



void
DSSC::parse_options(void)
{
  const ModelOptions& opts = SimulationInterface::get_options();

  _poisson_only = opts.get_option("poisson_only", _poisson_only);
}





void
DSSC::rebuild_equation_system(void)
{
  if (!_rebuild_eq_system) return;


  EquationSystems& equation_systems = get_equation_systems();


  // default is bcgsl
  ModelOptions& solveropts = get_solver_options();
  if (!solveropts.find_option("ksp_type"))
    solveropts["ksp_type"] = "bcgsl";

  // in 1D bcgs seems to work better than bcgsl
  const unsigned int dim = get_mesh().mesh_dimension();
  if ((dim == 1) && (solveropts["ksp_type"] == "bcgsl"))
    solveropts["ksp_type"] = "bcgs";


  if (solveropts.get_option("lin_abs_tol", -1.0) < 0)
    solveropts["lin_abs_tol"] = "1e-15";

  if (solveropts.get_option("nonlin_abs_tol", -1.0) < 0)
    solveropts["nonlin_abs_tol"] = "1e-15";


  // the coupled DD system
  TiberNonlinearSystem& system =
    *TiberNonlinearSystem::create(equation_systems,
      get_equation_system_name(), get_solver_options());
  
  system.attach_assembly_routine(assemble_system);


  system.add_variable("potential", libMeshEnums::FIRST);
  system.add_variable("fermi_n", libMeshEnums::FIRST);
  system.add_variable("fermi_I", libMeshEnums::FIRST);
  system.add_variable("fermi_I3", libMeshEnums::FIRST);
  system.add_variable("fermi_C", libMeshEnums::FIRST);


  // finally initialize the newly created system
  system.init();


  _rebuild_eq_system = false;

}




void
DSSC::do_init(void)
{

  _device = &get_environment().get_device();

  find_dirichlet_nodes();

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
      //ElectricalContact* contact = dynamic_cast<ElectricalContact*>(bd);
      //if (contact->is_real_contact())
      //{
      //  _boundary_currents[it->second] = 0.0;
      //  _voltages[it->second] = 0.0;
      //}
    }
  }
}




NumericVector<double>&
DSSC::solution_vector(void)
{
  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());
  
  return system.get_solution_vector();
}




void
DSSC::do_newton(void)
{

  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());


  system.set_options(get_solver_options());
  system.solve();
}





void
DSSC::find_dirichlet_nodes(void)
{
  SimulationEnvironment& env = get_environment();

  SimulationEnvironment::BoundaryNodeIterator it(env.boundary_nodes_begin());
  SimulationEnvironment::BoundaryNodeIterator end(env.boundary_nodes_end());
  
  for ( ; it != end; ++it)
  {
    ID id = it->second;

    Boundary* bd = env.get_boundary(id);

    DSSCContact* contact = NULL;
    if (bd != NULL)
      contact = dynamic_cast<DSSCContact*>(
          bd->get_boundary_properties(get_id()));

    if (contact != NULL)
    {
      //if ((contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
      //    || (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
      //    || (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET))
      //{
        _dirichlet_nodes[it->first] = bd;
      //}
    }
  }
}



/*
ID
DSSC::convert_variable_name_to_id(const string& variable_name) const
{
  ID id = INVALID_ID;
  
  // for an empty string we return immediately
  if (variable_name == "") return id;
  
  switch (variable_name[0])
  {
  case 'E':
      if (variable_name == "ElPotential")
        id = ELPOTENTIAL;
      break;
      
    case 'Q':
      if (variable_name == "QFermi_e")
        id = QFERMIE;
      break;

    default:
      break;
  }

  return id;
}

        


void
DSSC::get_solution_secure(const Elem* elem,
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

  DSSCModel* sc =
    dynamic_cast<DSSCModel*>(
        device.get_material(subdomain)->get_model(get_id()));

  assert(sc != NULL); 



  //sc->lock();
  sc->reinit(elem);

  fe->reinit(elem, &points);
    
  // Get the thermoelectric power
  sc->compute_thermoelectric_powers(); 
  double Pn_el =  sc->get_electron_thermoelectric_power();
  double Pp_el =  sc->get_hole_thermoelectric_power();
 
  vector<double> T_nodes = sc->get_temperature_at_nodes();

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
    double T = 0.0;
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

      T +=  phi[i][n] * T_nodes[i];    
 
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


    sc->set_coordinates(*(elem->get_node(n)));

    sc->set_potentials(u, en, ep);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_e(RealGradient(en_x, en_y, en_z));
    sc->set_grad_fermi_h(RealGradient(ep_x, ep_y, ep_z));

    sc->calculate_densities();

    sc->calculate_mobilities();
 
    sc->compute_thermoelectric_powers(); 

    double Pn =  sc->get_electron_thermoelectric_power();
    double Pp =  sc->get_hole_thermoelectric_power();
    
    double sigma_e = Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();

    double jnx = -sigma_e * (en_x + Pn_el * dT_x);
    double jny = -sigma_e * (en_y + Pn_el * dT_y);
    double jnz = -sigma_e * (en_z + Pn_el * dT_z);
    double jpx = -sigma_h * (ep_x + Pp_el * dT_x);
    double jpy = -sigma_h * (ep_y + Pp_el * dT_y);
    double jpz = -sigma_h * (ep_z + Pp_el * dT_z);

    

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
      values[n][SIGMAE] = sigma_e;

    if (ids.count(SIGMAH))
      values[n][SIGMAH] = sigma_h;

    if (ids.count(E))
      values[n][E] = e_field.size();

    if (ids.count(EX))
      values[n][EX] = e_field(0);

    if (ids.count(EY))
      values[n][EY] = e_field(1);

    if (ids.count(EZ))
      values[n][EZ] = e_field(2);

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
    
    if (ids.count(EJOULE))
      values[n][EJOULE] = ( jnx * jnx + jny * jny + jnz * jnz ) / sigma_e;

    if (ids.count(HJOULE))
      values[n][HJOULE] =  (jpx * jpx + jpy * jpy + jpz * jpz )/ sigma_h;
    
    if (ids.count(POWERNX))
      values[n][POWERNX] = (Pn * T + en) * jnx;
    
    if (ids.count(POWERNY))
      values[n][POWERNY] = (Pn * T + en) * jny;

    if (ids.count(POWERNZ))
      values[n][POWERNZ] = (Pn * T + en) * jnz;

    if (ids.count(POWERPX))
      values[n][POWERPX] = (Pp * T + ep) * jpx;
    
    if (ids.count(POWERPY))
      values[n][POWERPY] = (Pp * T + ep) * jpy;

    if (ids.count(POWERPZ))
      values[n][POWERPZ] = (Pp * T + ep) * jpz;
 
    if (ids.count(PN))
      values[n][PN] = Pn;
    
    if (ids.count(PP))
      values[n][PP] = Pp;

    if (ids.count(EPTSOURCE))
    {
      sc->compute_thermoelectric_power_gradient(); 
      RealGradient PnGrad =  sc->get_electron_thermoelectric_power_gradient();     
      values[n][EPTSOURCE] =  -T * ( PnGrad(0) * jnx + PnGrad(1) * jny + PnGrad(2) * jnz );
    }
    
    if (ids.count(HPTSOURCE))
    {
      sc->compute_thermoelectric_power_gradient(); 
      RealGradient PpGrad =  sc->get_hole_thermoelectric_power_gradient();
      values[n][HPTSOURCE] =  -T * ( PpGrad(0) * jpx + PpGrad(1) * jpy + PpGrad(2) * jpz );
    }

    if (ids.count(HRECOMB))
    { 
      vector<ID> rec_model_ids;
      int n_rec = sc->get_net_recombination_rate_IDs(rec_model_ids);
      double rec = 0.0;
      for (int i = 0; i < n_rec; i++)
	rec += sc->get_net_recombination_rate(rec_model_ids[i]);
      
      values[n][HRECOMB] = Constants::e * rec * (ep-en + T * (Pp-Pn));
    }
    

    set<ID>::iterator first(ids.begin());
    set<ID>::iterator it(ids.end());
    --it; 

    while (*it > MODELS )
    {
      if (*it > HEATMODELS)
      {
	double rec= sc->get_net_recombination_rate(*it - HEATMODELS);
	values[n][*it] = Constants::e * rec * (ep-en + T * (Pp-Pn));
      }
      else
	values[n][*it] = sc->get_net_recombination_rate(*it - MODELS);
      
      if (it == first)
        break;
      --it;

    }

   





  }

  //sc->unlock();
}

      


void
DSSC::get_solution_secure(const Elem* elem, const vector<Point>& p,
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

  DSSCModel* sc =
    dynamic_cast<DSSCModel*>(
        device.get_material(subdomain)->get_model(get_id()));
  assert(sc != NULL); 



  //sc->lock();
  sc->reinit(elem); //centroid

  fe->reinit(elem, &points);

    
  //Get the thermoelectric power
  sc->compute_thermoelectric_powers();
  double Pn_el =  sc->get_electron_thermoelectric_power();
  double Pp_el =  sc->get_hole_thermoelectric_power();
 
  vector<double> T_nodes = sc->get_temperature_at_nodes();

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
    double T = 0.0;
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

      T += phi[i][n] * T_nodes[i];

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
 

   
    sc->set_coordinates(p[n]);

    sc->set_potentials(u, en, ep);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_e(RealGradient(en_x, en_y, en_z));
    sc->set_grad_fermi_h(RealGradient(ep_x, ep_y, ep_z));


    sc->calculate_densities();

    sc->calculate_mobilities();
  
    sc->compute_thermoelectric_powers();

    double Pn =  sc->get_electron_thermoelectric_power();
    double Pp =  sc->get_hole_thermoelectric_power();

    double sigma_e = Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();


    double jnx = -sigma_e * (en_x + Pn_el * dT_x);
    double jny = -sigma_e * (en_y + Pn_el * dT_y);
    double jnz = -sigma_e * (en_z + Pn_el * dT_z);
    double jpx = -sigma_h * (ep_x + Pp_el * dT_x);
    double jpy = -sigma_h * (ep_y + Pp_el * dT_y);
    double jpz = -sigma_h * (ep_z + Pp_el * dT_z);

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
      values[n][SIGMAE] = sigma_e;

    if (ids.count(SIGMAH))
      values[n][SIGMAH] = sigma_h;

    if (ids.count(E))
      values[n][E] = e_field.size();

    if (ids.count(EX))
      values[n][EX] = e_field(0);

    if (ids.count(EY))
      values[n][EY] = e_field(1);

    if (ids.count(EZ))
      values[n][EZ] = e_field(2);
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

    if (ids.count(EJOULE))
      values[n][EJOULE] = ( jnx * jnx + jny * jny + jnz * jnz ) / sigma_e;
    
    if (ids.count(HJOULE))
      values[n][HJOULE] =  (jpx * jpx + jpy * jpy + jpz * jpz )/ sigma_h;
    
    if (ids.count(POWERNX))
      values[n][POWERNX] = (Pn * T + en) * jnx;
    
    if (ids.count(POWERNY))
      values[n][POWERNY] = (Pn * T + en) * jny;

    if (ids.count(POWERNZ))
      values[n][POWERNZ] = (Pn * T + en) * jnz;

    if (ids.count(POWERPX))
      values[n][POWERPX] = (Pp * T + ep) * jpx;
    
    if (ids.count(POWERPY))
      values[n][POWERPY] = (Pp * T + ep) * jpy;

    if (ids.count(POWERPZ))
      values[n][POWERPZ] = (Pp * T + ep) * jpz;

    if (ids.count(EPTSOURCE))
    {
      sc->compute_thermoelectric_power_gradient(); 
      RealGradient PnGrad =  sc->get_electron_thermoelectric_power_gradient();     
      values[n][EPTSOURCE] =  -T * ( PnGrad(0) * jnx + PnGrad(1) * jny + PnGrad(2) * jnz );
    }
    
    if (ids.count(HPTSOURCE))
    {
      sc->compute_thermoelectric_power_gradient(); 
      RealGradient PpGrad =  sc->get_hole_thermoelectric_power_gradient();
      values[n][HPTSOURCE] =  -T * ( PpGrad(0) * jpx + PpGrad(1) * jpy + PpGrad(2) * jpz );
    }

    if (ids.count(PN))
      values[n][PN] = Pn;
    
    if (ids.count(PP))
      values[n][PP] = Pp;

    if (ids.count(HRECOMB))
    { 
      vector<ID> rec_model_ids;
      int n_rec = sc->get_net_recombination_rate_IDs(rec_model_ids);
      double rec = 0.0;
      for (int i = 0; i < n_rec; i++)
	rec += sc->get_net_recombination_rate(rec_model_ids[i]);
          
      values[n][HRECOMB] = Constants::e * rec * (ep-en + T * (Pp-Pn));
    }
   

    set<ID>::iterator first(ids.begin());
    set<ID>::iterator it(ids.end());
    --it; 


    while (*it > MODELS )
    {
      if (*it > HEATMODELS)
      {
	double rec= sc->get_net_recombination_rate(*it - HEATMODELS);
	values[n][*it] = Constants::e * rec * (ep-en + T * (Pp-Pn));
      }
      else
	values[n][*it] = sc->get_net_recombination_rate(*it - MODELS);
      
      if (it == first)
        break;
      --it;

    }
   


  }

  //sc->unlock();
}
*/




void
DSSC::calculate_currents_rstf(void)
{
/*
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
  QGauss qrule(dim, libMeshEnums::FIFTH);
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

    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);

    
    fe->reinit(elem);

    sc->reinit(elem);
    
    //Get the temperature given the element
    //vector<double> T_nodes =  sc->get_temperature_at_nodes();

        
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
      sc->set_grad_fermi_e(phi0 * dEfn);
      sc->set_grad_fermi_h(phi0 * dEfp);

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
        {
          ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
                boundary->get_boundary_properties(get_id()));
          if (contact->is_real_contact())
            _boundary_currents[boundary] += j * dphi[n][qp];
        }

      }
    } // end loop over quadrature points
  } // end loop over elements
*/
}






void
DSSC::build_local_scaling(void)
{
/* 
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
  const double x0 = scaling.get_length_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
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

    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
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
      sc->set_electric_field(phi0 / x0 * e_field);
      sc->set_grad_fermi_e(RealGradient(0.0));
      sc->set_grad_fermi_h(RealGradient(0.0));

      sc->calculate_densities();
      sc->calculate_mobilities();

      double epsilon = sc->get_relative_permittivity();
      double l2_eps = JxW[qp] * l2 * epsilon;
      //double l2_eps = JxW[qp] * epsilon;

      double sigma_e = JxW[qp] *
        sc->get_electron_mobility() * sc->get_electron_density() / mu0;
      double sigma_h = JxW[qp] *
        sc->get_hole_mobility() * sc->get_hole_density() / mu0;


      double dn_dphi = sc->get_electron_density_derivative();
      double dp_dphi = sc->get_hole_density_derivative();
      double dNd_dphi = sc->get_ionized_donor_density_derivative();
      double dNa_dphi = sc->get_ionized_acceptor_density_derivative();

      double drho;
      drho = JxW[qp] * (dp_dphi - dNa_dphi - dn_dphi + dNd_dphi) * phi0 / C0;
      if (sc->is_dielectric())
        drho = 0.0;



      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //local_scaling_[elem->get_node(i)][0] += sigma_e * phi[i][qp];
        //local_scaling_[elem->get_node(i)][1] += sigma_h * phi[i][qp];
        local_scaling_[elem->get_node(i)][0] +=
          sigma_e * (dphi[i][qp] * dphi[i][qp]);
        local_scaling_[elem->get_node(i)][1] +=
          sigma_h * (dphi[i][qp] * dphi[i][qp]);
        local_scaling_[elem->get_node(i)][2] +=
          l2_eps * (dphi[i][qp] * dphi[i][qp]) -
          drho * phi[i][qp] * phi[i][qp];
      }


    } // end loop over quadrature points
  } // end loop over elements
*/
}



NumericVector<double>&
DSSC::do_get_solution_vector(void)
{
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  return system->get_solution_vector();
}



// implementation taken from libmesh equation_systems.C
void
DSSC::build_nodal_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{
 
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = get_solution_vector();

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

  int Pot = -1;
  if (variables.find("Potentials") != varend)
  {
    Pot = n_vars;
    n_vars += 5;
    legend.reserve(n_vars);
    legend.push_back("electric_potential");

    legend.push_back("el_chem_n");
    legend.push_back("el_chem_I");
    legend.push_back("el_chem_I3");
    legend.push_back("el_chem_C");
  }

  int Dens = -1;
  if (variables.find("Densities") != varend)
  {
    Dens = n_vars;
    n_vars += 4;
    legend.reserve(n_vars);
    legend.push_back("density_n");
    legend.push_back("density_I");
    legend.push_back("density_I3");
    legend.push_back("density_C");
  }

  int rho = -1;
  if (variables.find("TotalChargeDensity") != varend)
  {
    rho = n_vars;
    n_vars++;
    legend.reserve(n_vars);
    legend.push_back("total_charge_densitity");
  }

  int rec = -1;
  if (variables.find("Recombination") != varend)
  {
    rec = n_vars;
    n_vars++;
    legend.reserve(n_vars);
    legend.push_back("electron_recombination");
  }

  int gen = -1;
  if (variables.find("Generation") != varend)
  {
    gen = n_vars;
    n_vars++;
    legend.reserve(n_vars);
    legend.push_back("electron_generation");
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
  const unsigned int en_var = system->variable_number("fermi_n");
  const unsigned int eI_var = system->variable_number("fermi_I");
  const unsigned int eI3_var = system->variable_number("fermi_I3");
  const unsigned int eC_var = system->variable_number("fermi_C");
  
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;

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
      dof_map.dof_indices(elem, dof_indices_eI, eI_var);
      dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
      dof_map.dof_indices(elem, dof_indices_eC, eC_var);

      DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));
      assert(sc != NULL); 

      sc->reinit(elem);

      vector<Point> points(elem->n_nodes());
      for (int i = 0; i < elem->n_nodes(); i++)
        points[i] = elem->point(i);
      fe->reinit(elem, &points);

      assert(elem->n_nodes() == dof_indices_u.size());


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        double u  = phi0 * solution(dof_indices_u[n]);
        double en = phi0 * solution(dof_indices_en[n]);
        double eI = phi0 * solution(dof_indices_eI[n]);
        double eI3 = phi0 * solution(dof_indices_eI3[n]);
        double eC = phi0 * solution(dof_indices_eC[n]);

        RealGradient e_field(0);
        RealGradient grad_en(0);
        RealGradient grad_eI(0);
        RealGradient grad_eI3(0);
        RealGradient grad_eC(0);
        for (unsigned int i = 0; i < dof_indices_u.size(); i++)
        {
          e_field += dphi[i][n] * solution(dof_indices_u[n]);
          grad_en += dphi[i][n] * solution(dof_indices_en[n]);
          grad_eI += dphi[i][n] * solution(dof_indices_eI[n]);
          grad_eI3 += dphi[i][n] * solution(dof_indices_eI3[n]);
          grad_eC += dphi[i][n] * solution(dof_indices_eC[n]);
        }
        e_field *= -phi0;
        grad_en *= -phi0;
        grad_eI *= -phi0;
        grad_eI3 *= -phi0;
        grad_eC *= -phi0;



        // prepare for calculating local properties
        sc->set_coordinates(elem->point(n));


        sc->set_potentials(u, en, eI, eI3, eC);
        sc->set_electric_field(e_field);
        sc->set_grad_fermi_n(grad_en);
        sc->set_grad_fermi_I(grad_eI);
        sc->set_grad_fermi_I3(grad_eI3);
        sc->set_grad_fermi_C(grad_eC);

        sc->calculate_densities();
        sc->calculate_net_recombination_rate();


        assert (node_conn[elem->node(n)] != 0);
        double conn = static_cast<double>(node_conn[elem->node(n)]);

        unsigned int id = n_vars * elem->node(n);

        if (Dens != -1)
        {
          local[id + Dens] += sc->get_density_n() / conn;
          local[id + Dens + 1] += sc->get_density_I() / conn;
          local[id + Dens + 2] += sc->get_density_I3() / conn;
          local[id + Dens + 3] += sc->get_density_C() / conn;
        }

        if (rho != -1)
        {
          double nodal_val = sc->get_charge_density();
          local[id + rho] += nodal_val / conn;
        }

        if (rec != -1)
        {
          double nodal_val = sc->get_recombination_rate();
          local[id + rec] += nodal_val / conn;
        }

        if (gen != -1)
        {
          double nodal_val = sc->get_generation_rate();
          local[id + gen] += nodal_val / conn;
        }


        if (Pot != -1)
        {
          local[id + Pot] += u / conn;
          local[id + Pot + 1] += en / conn;
          local[id + Pot + 2] += eI / conn;
          local[id + Pot + 3] += eI3 / conn;
          local[id + Pot + 4] += eC / conn;
        }

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
DSSC::build_elemental_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{
/*
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

  int BandEdges = -1;
  if (variables.find("BandEdges") != varend)
  {
    const Elem* elem = *mesh.active_elements_begin();
    
    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(elem->subdomain_id())->get_model(get_id()));

    assert(sc != NULL); 

    sc->reinit(elem);
    sc->calculate_equilibrium_properties();

    const vector<double>& cb = sc->get_conduction_bands();
    const vector<double>& vb = sc->get_valence_bands();
    int n_bands = cb.size() + vb.size();
    legend.resize(legend.size() + n_bands);
    BandEdges = n_vars;

    for (unsigned int i = 0; i < cb.size(); i++)
    {
      ostringstream os;
      os << "CB" << i;
      legend[n_vars] = os.str();
      n_vars++;
    }

    for (unsigned int i = 0; i < vb.size(); i++)
    {
      ostringstream os;
      os << "VB" << i;
      legend[n_vars] = os.str();
      n_vars++;
    }
  }

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

  int GradFermiE = -1;
  if (variables.find("GradFermiE") != varend)
  {
    legend.resize(legend.size() + dim);
    GradFermiE = n_vars;
    switch (dim)
    {
      case 3:
        legend[GradFermiE + 2] = "grad_fermi_e_z";
        n_vars++;
      case 2:
        legend[GradFermiE + 1] = "grad_fermi_e_y";
        n_vars++;
        legend[GradFermiE + dim] = "mod_grad_fermi_e";
        n_vars++;
      default:
        legend[GradFermiE] = "grad_fermi_e_x";
        n_vars++;
    }
  }

  int GradFermiH = -1;
  if (variables.find("GradFermiH") != varend)
  {
    legend.resize(legend.size() + dim);
    GradFermiH = n_vars;
    switch (dim)
    {
      case 3:
        legend[GradFermiH + 2] = "grad_fermi_h_z";
        n_vars++;
      case 2:
        legend[GradFermiH + 1] = "grad_fermi_h_y";
        n_vars++;
        legend[GradFermiH + dim] = "mod_grad_fermi_h";
        n_vars++;
      default:
        legend[GradFermiH] = "grad_fermi_h_x";
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

  int dPn = -1;
  if (variables.find("GradPn") != varend)
  {
    legend.resize(legend.size() + dim);
    dPn = n_vars;
    switch (dim)
    {
      case 3:
        legend[dPn + 2] = "GradTepE_z";
        n_vars++;
      case 2:
        legend[dPn + 1] = "GradTepE_y";
        n_vars++;
        legend[dPn + dim] = "modGradTepE";
        n_vars++;
      default:
        legend[dPn] = "GradTepE_x";
        n_vars++;
    }
  }

  int dPp = -1;
  if (variables.find("GradPp") != varend)
  {
    legend.resize(legend.size() + dim);
    dPp = n_vars;
    switch (dim)
    {
      case 3:
        legend[dPp + 2] = "GradTepH_z";
        n_vars++;
      case 2:
        legend[dPp + 1] = "GradTepH_y";
        n_vars++;
        legend[dPp + dim] = "modGradTepH";
        n_vars++;
      default:
        legend[dPp] = "GradTepH_x";
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
  if (variables.find("PowerDensity") != varend)
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

    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL); 

    sc->reinit(elem);

    fe->reinit(elem);
    
  
    //Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_at_nodes();


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
    double oldu  = 0.0;
    double olden = 0.0;
    double oldep = 0.0;
    double T = 0.0;
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

      T += phi[i][0] * T_nodes[i];

      u  += phi[i][0] * solution(dof_indices_u[i]);
      en += phi[i][0] * solution(dof_indices_en[i]);
      ep += phi[i][0] * solution(dof_indices_ep[i]);
      oldu  += phi[i][0] * solution(dof_indices_u[i]);
      olden += phi[i][0] * solution(dof_indices_en[i]);
      oldep += phi[i][0] * solution(dof_indices_ep[i]);

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
    sc->set_old_potentials(phi0 * oldu, phi0 * olden, phi0 * oldep);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_e(RealGradient(en_x, en_y, en_z));
    sc->set_grad_fermi_h(RealGradient(ep_x, ep_y, ep_z));

    sc->calculate_densities();

    sc->calculate_mobilities();

    sc->compute_thermoelectric_powers();

    sc->compute_thermoelectric_power_gradient();

    double Pn = sc->get_electron_thermoelectric_power();
    double Pp = sc->get_hole_thermoelectric_power();
    RealGradient GradPn = sc->get_electron_thermoelectric_power_gradient();
    RealGradient GradPp = sc->get_hole_thermoelectric_power_gradient();

      
    double sigma_e = Constants::e * sc->get_electron_density() *
      sc->get_electron_mobility();
    double sigma_h = Constants::e * sc->get_hole_density() *
      sc->get_hole_mobility();

    unsigned int id = n_vars * elem_number;

    if (BandEdges != -1)
    {
      sc->calculate_equilibrium_properties();
      const vector<double>& cb = sc->get_conduction_bands();
      const vector<double>& vb = sc->get_valence_bands();

      for (unsigned int i = 0; i < cb.size(); i++)
        results[id + BandEdges + i] = cb[i];

      for (unsigned int i = 0; i < vb.size(); i++)
        results[id + BandEdges + cb.size() + i] = vb[i];
    }


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


    if (GradFermiE != -1)
    {
      switch (dim)
      {
        case 3:
          results[id + GradFermiE + 2] = en_z;
        case 2:
          results[id + GradFermiE + 1] = en_y;
          results[id + GradFermiE + dim] =  sqrt(en_x * en_x + en_y * en_y + en_z * en_z);
        default:
          results[id + GradFermiE] = en_x;
      }
    }


    if (GradFermiH != -1)
    {
      switch (dim)
      {
        case 3:
          results[id + GradFermiH + 2] = ep_z;
        case 2:
          results[id + GradFermiH + 1] = ep_y;
          results[id + GradFermiH + dim] =
            sqrt(ep_x * ep_x + ep_y * ep_y + ep_z * ep_z);
        default:
          results[id + GradFermiH] = ep_x;
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
      double jx = -sigma_e * (en_x + Pn * dT_x);
      double jy = -sigma_e * (en_y + Pn * dT_y);
      double jz = -sigma_e * (en_z + Pn * dT_z);
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
      double jx = -sigma_h * (ep_x + Pp * dT_x);
      double jy = -sigma_h * (ep_y + Pp * dT_y);
      double jz = -sigma_h * (ep_z + Pp * dT_z);
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
      double jx = -sigma_e * (en_x + Pn * dT_x) - sigma_h * (ep_x + Pp * dT_x);
      double jy = -sigma_e * (en_y + Pn * dT_y) - sigma_h * (ep_y + Pp * dT_y);
      double jz = -sigma_e * (en_z + Pn * dT_z) - sigma_h * (ep_z + Pp * dT_z);
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


    if (dPn != -1)
    {
      double dPn_x = GradPn(0);
      double dPn_y = GradPn(1);
      double dPn_z = GradPn(2);
      switch (dim)
      {
        case 3:
          results[id + dPn + 2] = dPn_z;
        case 2:
          results[id + dPn + 1] = dPn_y;
          results[id + dPn + dim] = sqrt(dPn_x * dPn_x + dPn_y * dPn_y + dPn_z * dPn_z);
        default:
          results[id + dPn] = dPn_x;
      }
    }

    if (dPp != -1)
    {
      double dPp_x = GradPp(0);
      double dPp_y = GradPp(1);
      double dPp_z = GradPp(2);
      switch (dim)
      {
        case 3:
          results[id + dPp + 2] = dPp_z;
        case 2:
          results[id + dPp + 1] = dPp_y;
          results[id + dPp + dim] = sqrt(dPp_x * dPp_x + dPp_y * dPp_y + dPp_z * dPp_z);
        default:
          results[id + dPp] = dPp_x;
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


  


*/
}



void
DSSC::build_integrated_quantities(const set<string>& names,
    vector<double>& values)
{
/*
  const set<string>::const_iterator varend(names.end());

  if ((names.find("ContactCurrents") != varend) ||
      (names.find("current") != varend))
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
*/
}


void
DSSC::calculate_currents(void)
{
  calculate_currents_rstf();
}


void
DSSC::build_integrated_quantities_description(
    const std::set<std::string>& names,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  const set<string>::const_iterator varend(names.end());

  if ((names.find("ContactCurrents") != varend) ||
      (names.find("current") != varend))
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
DSSC::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_potential_scaling();
}






void
DSSC::assemble_system(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{
  _this->do_assembly(x, residual, jacobian);
}






void
DSSC::do_assembly(const NumericVector<Number>& x,
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


  BoundaryNodeList& dirichlet_nodes = _dirichlet_nodes;


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
  // scaling for electrons
  double C0_e = C0;
  // scaling for I
  double C0_I = C0;
  // scaling for I3
  double C0_I3 = C0;
  // scaling for C
  double C0_C = C0;
 
  if (do_local_scaling_)
    C0_e = C0_I = C0_I3 = C0_C = 1.0;


  // scaling for recombination rates
  double R0_e = C0_e / scaling.get_time_scaling();
  double R0_I = C0_I / scaling.get_time_scaling();
  double R0_I3 = C0_I3 / scaling.get_time_scaling();
  //double R0_C = C0_C / scaling.get_time_scaling();


  const DofMap& dof_map = system.get_dof_map();
  
  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_n");
  const unsigned int eI_var = system.variable_number("fermi_I");
  const unsigned int eI3_var = system.variable_number("fermi_I3");
  const unsigned int eC_var = system.variable_number("fermi_C");
  
  FEType fe_type = system.variable_type(u_var);

  libMeshEnums::Order integration_order = libMeshEnums::FIFTH;

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

  DenseSubMatrix<Number>
    Kuu(Ke), Kun(Ke), KuI(Ke), KuI3(Ke), KuC(Ke),
    Knu(Ke), Knn(Ke), KnI(Ke), KnI3(Ke),  // KnC = 0
    KIu(Ke), KIn(Ke), KII(Ke), KII3(Ke),  // KIC = 0
    KI3u(Ke),KI3n(Ke),KI3I(Ke),KI3I3(Ke), // KI3C = 0
    KCu(Ke), /* =0 */ /* =0 */ /* =0 */  KCC(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fn(Fe),
    FI(Fe),
    FI3(Fe),
    FC(Fe);

  DenseSubVector<Number>
    Xu(X),
    Xn(X),
    XI(X),
    XI3(X),
    XC(X);


  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;

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
    dof_map.dof_indices(elem, dof_indices_eI, eI_var);
    dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
    dof_map.dof_indices(elem, dof_indices_eC, eC_var);

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
    //         -                       -          -  -
    //        | Kuu  Kun  KuI  KuI3 KuC |        | Fu |
    //        | Knu  Knn  KnI  KnI3  0  |        | Fn |
    //   Ke = | KIu  KIn  KII  KII3  0  |;  Fe = | FI |
    //        | KI3u KI3n KI3I KI3I3 0  |        | FI3|
    //        | KCu   0    0     0  KCC |        | FC |
    //         -                       -          -  -
    //
    Kuu.reposition(0, 0, n_dofs, n_dofs);
    Kun.reposition(0, n_dofs, n_dofs, n_dofs);
    KuI.reposition(0, 2 * n_dofs, n_dofs, n_dofs);
    KuI3.reposition(0, 3 * n_dofs, n_dofs, n_dofs);
    KuC.reposition(0, 4 * n_dofs, n_dofs, n_dofs);
    //
    Knu.reposition(n_dofs, 0, n_dofs, n_dofs);
    Knn.reposition(n_dofs, n_dofs, n_dofs, n_dofs);
    KnI.reposition(n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KnI3.reposition(n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    //
    KIu.reposition(2 * n_dofs, 0, n_dofs, n_dofs);
    KIn.reposition(2 * n_dofs, n_dofs, n_dofs, n_dofs);
    KII.reposition(2 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KII3.reposition(2 * n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    //
    KI3u.reposition(3 * n_dofs, 0, n_dofs, n_dofs);
    KI3n.reposition(3 * n_dofs, n_dofs, n_dofs, n_dofs);
    KI3I.reposition(3 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KI3I3.reposition(3 * n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    //
    KCu.reposition(4 * n_dofs, 0, n_dofs, n_dofs);
    KCC.reposition(4 * n_dofs, 4 * n_dofs, n_dofs, n_dofs);
    //
    Fu.reposition(0, n_dofs);
    Fn.reposition(n_dofs, n_dofs);
    FI.reposition(2 * n_dofs, n_dofs);
    FI3.reposition(3 * n_dofs, n_dofs);
    FC.reposition(4 * n_dofs, n_dofs);
    //
    Xu.reposition(0, n_dofs);
    Xn.reposition(n_dofs, n_dofs);
    XI.reposition(2 * n_dofs, n_dofs);
    XI3.reposition(3 * n_dofs, n_dofs);
    XC.reposition(4 * n_dofs, n_dofs);



    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);
    sc->reinit(elem);
    

    // Get the temperature given the element
    //vector<double> T_nodes = sc->get_temperature_at_nodes();
   

    vector<vector<double> > local_scaling(elem->n_nodes(), vector<double>(3, 1));
    if (do_local_scaling_)
    {
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        local_scaling[n] = local_scaling_[elem->get_node(n)];
      }
    }


    bool is_electrolyte = sc->is_electrolyte();
    bool is_TiO2 = sc->is_TiO2();



    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real u  = 0.0;
      Real en = 0.0;
      Real eI = 0.0;
      Real eI3 = 0.0;
      Real eC = 0.0;
      RealGradient e_field(0);
      RealGradient grad_en(0);
      RealGradient grad_eI(0);
      RealGradient grad_eI3(0);
      RealGradient grad_eC(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * Xu(i);
        en += phi[i][qp] * Xn(i);
        eI += phi[i][qp] * XI(i);
        eI3 += phi[i][qp] * XI3(i);
        eC += phi[i][qp] * XC(i);
        e_field += dphi[i][qp] * Xu(i);
        grad_en += dphi[i][qp] * Xn(i);
        grad_eI += dphi[i][qp] * XI(i);
        grad_eI3 += dphi[i][qp] * XI3(i);
        grad_eC += dphi[i][qp] * XC(i);
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);

      double grad_fac = phi0 / x0;
      sc->set_electric_field(grad_fac * e_field);
      sc->set_grad_fermi_n(grad_fac * grad_en);
      sc->set_grad_fermi_I(grad_fac * grad_eI);
      sc->set_grad_fermi_I3(grad_fac * grad_eI3);
      sc->set_grad_fermi_C(grad_fac * grad_eC);

      // calculate all local properties
      sc->calculate_densities();
      sc->calculate_net_recombination_rate();

      double n_e = sc->get_density_n();
      double n_I = sc->get_density_I();
      double n_I3 = sc->get_density_I3();
      double n_C = sc->get_density_C();
  

    
      double epsilon = sc->get_relative_permittivity();
      double l2_eps = l2 * epsilon;

      double R = 0.0;
      if ((is_TiO2 && is_electrolyte) && !poisson_only())
        R = sc->get_net_recombination_rate();
      

      double mu_n = sc->get_mobility_n();
      double mu_I = sc->get_mobility_I();
      double mu_I3 = sc->get_mobility_I3();
      double mu_C = sc->get_mobility_C();


      // the jacobian x weight x scaling
      double J = JxW[qp];


      // NOTE: sigma_e = mu_e * n is the electron conductivity
      double sigma_n = mu_n * n_e / (mu0 * C0_e);
      double sigma_I = mu_I * n_I / (mu0 * C0_I);
      double sigma_I3 = mu_I3 * n_I3 / (mu0 * C0_I3);
      double sigma_C = mu_C * n_C / (mu0 * C0_C);

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
          
          Kuu(i,j) += l2_eps * laplace / local_scaling[i][2];
          
          if (!poisson_only())
          {
            if (is_TiO2)
              Knn(i,j) += sigma_n * laplace / local_scaling[i][0];

            if (is_electrolyte)
            {
              KII(i,j) += sigma_I * laplace / local_scaling[i][1];

              KI3I3(i,j) += sigma_I3 * laplace / local_scaling[i][1];

              KCC(i,j) += sigma_C * laplace / local_scaling[i][1];
            }
          }
        }

        if (!is_TiO2 || poisson_only())
          Knn(i,i) += 1;
        
        if (!is_electrolyte || poisson_only())
        {
          KII(i,i) += 1;
          KI3I3(i,i) += 1;
          KCC(i,i) += 1;
        }
      }
      
      // 
      // for jacobian compute the other contributions
      // 
      if (jacobian != NULL)
      {
        double dn_dphi = sc->get_density_derivative_n();
        double dI_dphi = sc->get_density_derivative_I();
        double dI3_dphi = sc->get_density_derivative_I3();
        double dC_dphi = sc->get_density_derivative_C();

        double drho = dC_dphi - (dn_dphi + dI_dphi + dI3_dphi);
        drho *= phi0 / C0;


        double dR[4];

        if ((R == 0.0) || poisson_only() || !(is_TiO2 && is_electrolyte))
          dR[0] = dR[1] = dR[2] = dR[3] = 0.0;
        else
        {
          double dR_dn = sc->get_net_recombination_rate_derivatives()[0];
          double dR_dI = sc->get_net_recombination_rate_derivatives()[1];
          double dR_dI3 = sc->get_net_recombination_rate_derivatives()[2];

          dR[1] = -dR_dn * dn_dphi * phi0;
          dR[2] = -dR_dI * dI_dphi * phi0;
          dR[3] = -dR_dI3 * dI3_dphi * phi0;
          dR[0] = -(dR[1] + dR[2] + dR[3]);
        }

        // d(sigma_n)/du * element-jacobian
        // sigma_n = mu_n * n means the conductivity of electrons
        // the factor phi_0 comes from the derivative with respect to the potential
        double dsigma_n = J * phi0 / (mu0 * C0_e) * mu_n * dn_dphi;
        double dsigma_I = J * phi0 / (mu0 * C0_I) * mu_I * dI_dphi;
        double dsigma_I3 = J * phi0 / (mu0 * C0_I3) * mu_I3 * dI3_dphi;
        double dsigma_C = J * phi0 / (mu0 * C0_C) * mu_C * dC_dphi;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double lap_e = (dphi[i][qp] * grad_en) / local_scaling[i][0];
          double lap_I = (dphi[i][qp] * grad_eI) / local_scaling[i][1];
          double lap_I3 = (dphi[i][qp] * grad_eI3) / local_scaling[i][1];
          double lap_C = (dphi[i][qp] * grad_eC) / local_scaling[i][1];
          double dsigma_n_x_lap = dsigma_n * lap_e;
          double dsigma_I_x_lap = dsigma_I * lap_I;
          double dsigma_I3_x_lap = dsigma_I3 * lap_I3;
          double dsigma_C_x_lap = dsigma_C * lap_C;

          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            //
            // NOTE: we do not make a loop over l, but use precalculated
            // gradients

            double dsigma_n_x_phi = dsigma_n_x_lap * phi[j][qp];
            double dsigma_I_x_phi = dsigma_I_x_lap * phi[j][qp];
            double dsigma_I3_x_phi = dsigma_I3_x_lap * phi[j][qp];
            double dsigma_C_x_phi = dsigma_C_x_lap * phi[j][qp];

            if (is_TiO2 && !poisson_only())
            {
              Knu(i,j) += dsigma_n_x_phi;
              Knn(i,j) -= dsigma_n_x_phi;
            }

            if (is_electrolyte && !poisson_only())
            {
              KIu(i,j) += dsigma_I_x_phi;
              KII(i,j) -= dsigma_I_x_phi;

              KI3u(i,j) += dsigma_I3_x_phi;
              KI3I3(i,j) -= dsigma_I3_x_phi;

              KCu(i,j) += dsigma_C_x_phi;
              KCC(i,j) -= dsigma_C_x_phi;
            }




            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            Kuu(i,j) -= drho * phi_i_x_phi_j / local_scaling[i][2];

            if (is_TiO2)
            {
              Kun(i,j) -= phi0 * dn_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_e;

              if (!poisson_only())
              {
                Knu(i,j) -= dR[0] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                Knn(i,j) -= dR[1] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
              }
            }

            if (is_electrolyte)
            {
              KuI(i,j) -= phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I;
              KuI3(i,j) -= phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I3;
              KuC(i,j) -= phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_C;

              if (!poisson_only())
              {
                if (is_TiO2)
                {
                  KnI(i,j) -= dR[2] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                  KnI3(i,j) -= dR[3] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;

                  KIn(i,j) -= 1.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                  KI3n(i,j) += 0.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][1] / R0_I3;
                }

                KIu(i,j) -= 1.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII(i,j) -= 1.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII3(i,j) -= 1.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;

                KI3u(i,j) += 0.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][1] / R0_I3;
                KI3I(i,j) += 0.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][1] / R0_I3;
                KI3I3(i,j) += 0.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][1] / R0_I3;
              }
            }
          }

        }

      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        double J_x_rho = J * sc->get_charge_density() / C0;

        // net recombination rate
        double J_x_R = 0.0;
        if ((is_TiO2 && is_electrolyte) && !poisson_only())
          J_x_R = J * R;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double net_recomb = J_x_R * phi[i][qp] / local_scaling[i][0];
          
          Fu(i) -= J_x_rho * phi[i][qp] / local_scaling[i][2];
          
          Fn(i) -= net_recomb / R0_e;

          // TODO factors
          FI(i) -= 1.5 * net_recomb / R0_I;
          FI3(i) += 0.5 * net_recomb / R0_I3;
        }
      }

    } // end loop over quadrature points
   


    set<unsigned int> nodes_on_boundary_sides;
 
    //cerr << Knn << endl;
    
///*
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

        DSSCContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<DSSCContact*>(
              boundary->get_boundary_properties(get_id()));


        fe_face->reinit(elem, s);

/*
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

        }
*/


/*
        // for von Neumann or mixed type boundary conditions
        vector<double> coeff(3, 0.0);
        vector<double> value(3, 0.0);

        // the derivatives
        vector<vector<double> > dcoeff(3, vector<double>(3, 0.0));
        vector<vector<double> > dvalue(3, vector<double>(3, 0.0));

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
              sc->set_grad_fermi_e(phi0 / x0 * grad_en);
              sc->set_grad_fermi_h(phi0 / x0 * grad_ep);
              sc->calculate_densities();
              sc->calculate_mobilities();

              // we put the phi0 here for convenience
              double sigma_e = phi0 * sc->get_electron_mobility() *
                sc->get_electron_density();
              double sigma_h = phi0 * sc->get_hole_mobility() *
                sc->get_hole_density();

              sc->compute_thermoelectric_powers();
              double Pn =  sc->get_electron_thermoelectric_power();
              double Pp =  sc->get_hole_thermoelectric_power();

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
              double scale_u = J * phi0 / x0 / C0;
              double scale_n = J * phi0 / x0 / mu0;

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                double fac_u = scale_u / local_scaling[i][2];

                for (unsigned int j = 0; j < n_dofs; j++)
                {

                  Real phi_i_x_phi_j =
                    phi_face[i][qp] * phi_face[j][qp];

                  if (coupling & POISSON)
                  {
                    Kuu(i,j) -= fac_u * dvalue[0][0] * phi_i_x_phi_j;

                    //if (coupling & ECURRENT)
                    //  Kun(i,j) -= fac_u * dvalue[0][1] * phi_i_x_phi_j;

                    //if (coupling & HCURRENT)
                    //  Kup(i,j) -= fac_u * dvalue[0][2] * phi_i_x_phi_j;
                  }

                  //if (coupling & ECURRENT)
                  //  Knn(i,j) += coeff[1] * phi_i_x_phi_j;

                  //if (coupling & HCURRENT)
                  //  Kpp(i,j) += coeff[2] * phi_i_x_phi_j;
                }
              }
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
*/
          nodes_on_boundary_sides.insert(elem->node(s));
/*
          // s is the node of the element lying on the boundary
          Real u  = Xu(s);
          Real en = Xn(s);
          Real eI = XI(s);
          Real eI3 = XI3(s);
          Real eC = XC(s);

          // calculate densities etc.
          sc->set_coordinates(elem->point(s));
          sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);
          
          RealGradient e_field(0.0);
          double grad_en = 0.0;
          double grad_eI = 0.0;
          double grad_eI3 = 0.0;
          double grad_eC = 0.0;
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            grad_en += dphi_face[n][0](0)* Xn(n);
            grad_eI += dphi_face[n][0](0) * XI(n);
            grad_eI3 += dphi_face[n][0](0) * XI3(n);
            grad_eC += dphi_face[n][0](0) * XC(n);
            e_field(0) += dphi_face[n][0](0) * Xu(n);

          }
          
          sc->set_electric_field(phi0 / x0 * e_field);
          sc->set_grad_fermi_n(phi0 / x0 * RealGradient(grad_en, 0.0, 0.0));
          sc->set_grad_fermi_I(phi0 / x0 * RealGradient(grad_eI, 0.0, 0.0));
          sc->set_grad_fermi_I3(phi0 / x0 * RealGradient(grad_eI3, 0.0, 0.0));
          sc->set_grad_fermi_C(phi0 / x0 * RealGradient(grad_eC, 0.0, 0.0));
          sc->calculate_densities();

          // we put the phi0 here for convenience
          double sigma_e = phi0 * sc->get_electron_mobility() *
            sc->get_electron_density();
          double sigma_h = phi0 * sc->get_hole_mobility() *
            sc->get_hole_density();

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
            double val_uu = dvalue[0][0] / x0 / C0 * phi0;
            double val_nn = dvalue[1][1] / mu0 / C0 * phi0;
            double val_pp = dvalue[2][2] / mu0 / C0 * phi0;

            if (coupling & POISSON)
              Kuu(s,s) -= val_uu / local_scaling[s][2];

            if (coupling & ECURRENT)
              Knn(s,s) -= val_nn / local_scaling[s][0];

            if (coupling & HCURRENT)
              Kpp(s,s) -= val_pp / local_scaling[s][1];
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
*/
      }
    } // end loop over element sides
//*/


    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    // NOTE: this changes dof_indices that's why the application of
    //       Dirichlet type BCs needs special care
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

///*
    //
    // now as last thing we apply Dirichlet type Bcs
    //
    BoundaryNodeList::const_iterator node_it;
    const BoundaryNodeList::const_iterator end =
      dirichlet_nodes.end();
//    if (Ke.m() == n_dofs_tot)
    {
      // no constrained nodes, so everything is easy

      // loop over all nodes and check if it is a dirichlet type node
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        node_it = dirichlet_nodes.find(elem->get_node(i));
        if (node_it != end)
        {
          Boundary* bd = node_it->second;
          DSSCContact* contact = dynamic_cast<DSSCContact*>(
              bd->get_boundary_properties(get_id()));
          //contact->set_material(sc);
          //contact->set_normal_fluxes(
          //    nodal_flux_n[elem->node(i)], nodal_flux_p[elem->node(i)]);

          // we only impose Dirichlet type BCs if the node has an associated
          // boundary side
          if (nodes_on_boundary_sides.find(elem->node(i)) !=
              nodes_on_boundary_sides.end())
          {
            double val = contact->get_potential() / phi0;
            Ke.condense(i + n_dofs, i + n_dofs, -val, Fe);
          }
/*
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
*/
        }
      }
    }
/*
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
*/


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
  {
    jacobian->close();
    jacobian->print_matlab("J.m");
  }
  else
  {
    residual->close();
    residual->print_matlab("F.m");
  }

  
  perf_log.stop_event("assembly");

} 



