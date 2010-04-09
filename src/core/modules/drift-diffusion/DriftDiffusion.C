// $Id$

// module includes
#include "DriftDiffusion.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "Scaling.h"
#include "Material.h"
#include "Boundary.h"
#include "ElectricalContact.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "TiberNonlinearSystem.h"
#include "SolveFailedException.h"
#include "Variable.h"

// TODO should be replaced by boost methods
#include "gzstream.h"

// libmesh includes
#include "node.h"
#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "equation_systems.h"
#include "mesh_refinement.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

#include "DataOutput.h"
#include "Messages.h"

// C++ includes
#include <fstream>

//
// Module interface
//

TIBER_MODULE(DriftDiffusion, MODULE_NAME)

namespace
{
  int __private_counter;
}


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
    quadrature_type(QTRAP),
    integration_order(libMeshEnums::FIFTH),
    solver_method(NEWTON),
    max_gummel_iterations(5),
    scaling_type(Scaling::UNITS),
    coupling(FULLYCOUPLED),
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
    quadrature_type(rhs.quadrature_type),
    integration_order(rhs.integration_order),
    solver_method(rhs.solver_method),
    max_gummel_iterations(rhs.max_gummel_iterations),
    scaling_type(rhs.scaling_type),
    coupling(rhs.coupling),
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
    quadrature_type = rhs.quadrature_type;
    integration_order = rhs.integration_order;
    solver_method = rhs.solver_method;
    max_gummel_iterations = rhs.max_gummel_iterations;
    scaling_type = rhs.scaling_type;
    coupling = rhs.coupling;
    current_calculation = rhs.current_calculation;
    exact_newton = rhs.exact_newton;
  }
  return *this;
}





DriftDiffusion::DriftDiffusion(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true),
    _useparticle('b')
{
}




DriftDiffusion::~DriftDiffusion(void)
{
  cleanup_solver();
}




PhysicalModel*
DriftDiffusion::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  string modelname;

  modelname = options.get_option("model", "default");

  if ((modelname == "unstrained") || (modelname == "strained"))
  {
    ostringstream os;
    os << "drift-diffusion model '" << modelname << "' is deprecated."
      << endl << "Use 'model = default' instead or don't specify explicitly."
      << endl;
    Messages::warning(os.str());
    modelname = "default";
  }

  DriftDiffusionProperties* model =
    DriftDiffusionProperties::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "DriftDiffusion: No such physical model: " + modelname);

  // we need non-const pointer to DriftDiffusion
  model->set_driftdiffusion(const_cast<DriftDiffusion*>(this));

  return model;
}



PhysicalModel*
DriftDiffusion::create_boundary_model(const ModelOptions& options,
    const Material* material_A, const Material* material_B) const
{

  return NULL;
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

  // we calculate in cm!
  double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  if (type == Scaling::NONE)
  {
    get_scaling().set_scaling_type(type);
    get_scaling().set_potential_scaling(1);
    get_scaling().set_length_scaling(1);
    get_scaling().set_mobility_scaling(1);
    get_scaling().set_density_scaling(1);

    // We don't have to do anything in this case
    return;
  }


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
          get_physical_model(elem->subdomain_id()));

    sc->set_coordinates(elem->centroid());
    sc->set_potentials(sc->get_equilibrium_fermi_level());
    sc->set_electric_field(RealGradient(0));
    sc->set_grad_fermi_e(RealGradient(0));
    sc->set_grad_fermi_h(RealGradient(0));
    sc->reinit(elem);

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
  double level = -Ef_n / phi0;

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

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
  double level = -Ef_p / phi0;

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

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

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

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

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_e");
  const unsigned int ep_var = system->variable_number("fermi_h");

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;


  MeshBase& mesh = get_mesh();
  unsigned int dim = mesh.mesh_dimension();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    for (size_t s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(elem->top_parent(), s);

      if (env.is_boundary(side))
      {
        Boundary* boundary = env.get_boundary(side);
        if (boundary == NULL)
          continue;

        ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
            boundary->get_boundary_properties(get_id()));

        if (contact != NULL)
        {
          if ((contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
              || (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
              || (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET))
          {
            AutoPtr<Elem> side_el(elem->build_side(s));
            if (dim == 1)
            {
              dof_indices_u = vector<unsigned int>(1,
                  elem->get_node(s)->dof_number(system->number(), u_var, 0));

              dof_indices_en = vector<unsigned int>(1,
                  elem->get_node(s)->dof_number(system->number(), en_var, 0));

              dof_indices_ep = vector<unsigned int>(1,
                  elem->get_node(s)->dof_number(system->number(), ep_var, 0));
            }
            else
            {
              dof_map.dof_indices(side_el.get(), dof_indices_u, u_var);
              dof_map.dof_indices(side_el.get(), dof_indices_en, en_var);
              dof_map.dof_indices(side_el.get(), dof_indices_ep, ep_var);
            }

            for (unsigned int i = 0; i < side_el->n_nodes(); i++)
            {
              if (contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
                _dirichlet_dofs.insert(dof_indices_u[i]);

              if (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
                _dirichlet_dofs.insert(dof_indices_en[i]);

              if (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET)
                _dirichlet_dofs.insert(dof_indices_ep[i]);
            }
          }
        }
      }
    }
  }



  /*{
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
  }*/
}



void
DriftDiffusion::find_dielectric_boundary_nodes(void)
{
  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* el = *it;

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          get_physical_model(el->subdomain_id()));

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
                get_physical_model(el->neighbor(s)->subdomain_id()));


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
    //get_equation_systems().delete_system(get_equation_system_name());
    _rebuild_eq_system = true;
  }
}




void
DriftDiffusion::cleanup_solver(void)
{

  // erase boundary current data structure
  _boundary_currents.erase(_boundary_currents.begin(),
      _boundary_currents.end());

  reset_solver();
}




void
DriftDiffusion::do_solve(void)
{
  __private_counter = 0;

  string filename = get_option("load_state", "");
  if (filename != "")
  {
    compute_scaling(get_my_options().scaling_type);
    Messages::info(get_name() + ": Loading state from " + filename);
    load_data(filename);
    get_options().set_option("load_state", "");
    if (!get_option("solve_after_load", false)) return;
  }

  // rebuild the system if needed
  //rebuild_equation_system();

  // set a static pointer to ourselves
  // this is needed in the static assembly routine
  _this = this;

  parse_options();


  bool equilibrium = true;
  {
    ContactData::iterator it(_voltages.begin());
    const ContactData::iterator end(_voltages.end());
    for ( ; it != end; ++it)
    {
      ElectricalContact* bd = static_cast<ElectricalContact*>(
          it->first->get_boundary_properties(get_id()));

      double voltage = bd->get_simulation_voltage();

      _voltages[it->first] = voltage;

      if (voltage != 0.0)
        equilibrium = false;
    }
  }


  if (!equilibrium_done())
  {
    solve_equilibrium();

    bool save = get_option("save_state", false);
    if (save)
    {
      string file = get_output_directory() + "/" +
      get_output_filename_prefix() + "_equilibrium.tsv";
      save_data(file);
    }

    if (do_local_scaling_)
      build_local_scaling();

    // if we would repeat the equilibrium simulation, we can stop now
    if (equilibrium)
      return;
  }

  // set the old solution
  EquationSystems& es = get_equation_systems();
  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());
  get_solution_vector().close();
  system.get_vector("old_sol") = get_solution_vector();

  int coupling = get_my_options().coupling;




  //set_dirichlet_bc();

  try
  {
    switch (_options.solver_method)
    {
      case GUMMEL:
        //solve_gummel();
        break;
      default: // Newton method
        do_newton();
        break;
    }
  }
  catch (SolverException& e)
  {
    string msg = "solve failed (" +
        string(e.what()) + ")";
    throw SolveFailedException(msg);
  }

  // NOTE we calculate the local scaling factors AFTER, because otherwise
  // new results of other coupled models (thermal, quantum) may disturb
  // the calculation
  if (do_local_scaling_)
    build_local_scaling();

  get_my_options().coupling = coupling;

  // calculate the currents to print them on screen
  calculate_currents();


  bool field_emission = false;

  ContactData::iterator it(_boundary_currents.begin());
  const ContactData::iterator end(_boundary_currents.end());
  for ( ; it != end; ++it)
  {
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(
          it->first->get_boundary_properties(get_id()));
    if (cnt->has_field_emission())
    {
      field_emission = true;
      break;
    }
  }

  int width = 20;
  {
    Messages::newline();
    ostringstream os;
    os << "contact name:";
    os.width(width - os.tellp());
    os << "";
    os << "contact voltage:";
    os.width(2 * width - os.tellp());
    os << "";
    os << "contact current:";
    if (field_emission)
    {
      os.width(3 * width - os.tellp());
      os << "";
      os << "field emission current:";
    }
    os << endl;
    Messages::info(os.str());
  }

  for (it = _boundary_currents.begin(); it != end; ++it)
  {
    ostringstream os;
    os << setprecision(6);
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(it->first->get_boundary_properties(get_id()));
    os << it->first->get_name();
    if (cnt->get_simulation_voltage() < 0)
      os.width(width - os.tellp() - 1);
    else
      os.width(width - os.tellp());
    os << "";
    os << cnt->get_simulation_voltage();

    if (it->second < 0)
      os.width(2 * width - os.tellp() - 1);
    else
      os.width(2 * width - os.tellp());
    os << "";
    os << it->second * it->first->get_area_factor();

    if (cnt->has_field_emission())
    {
      double iem = cnt->get_field_emission_current();
      if (iem < 0)
        os.width(3 * width - os.tellp() - 1);
      else
        os.width(3 * width - os.tellp());
      os << "";
      os << iem * it->first->get_area_factor();
    }
    Messages::info(os.str());
  }

  bool save = get_option("save_state", false);
  if (save)
  {
    string file = get_output_directory() + "/" +
      get_output_filename() + ".tsv";
    save_data(file);
  }
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
  compute_scaling(get_my_options().scaling_type);


  ModelOptions& solveropts = get_solver_options();
  int max_it = solveropts.get_option("nonlin_max_it", 15);
  solveropts.set_option("nonlin_max_it", 150);

  int coupling = get_my_options().coupling;
  get_my_options().coupling = POISSON;

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


  if (do_local_scaling_)
    build_local_scaling();

  const ModelOptions& opts = get_options();
  if (opts.find_option("el_qfermi_level"))
    set_electron_fermi_level(opts.get_option("el_qfermi_level", 0.0));

  if (opts.find_option("hl_qfermi_level"))
    set_hole_fermi_level(opts.get_option("hl_qfermi_level", 0.0));



  try
  {
    Messages::info("Solving equilibrium");

    do_newton();

    Messages::info("Equilibrium done");
  }
  catch (runtime_error& e)
  {
    ostringstream os;
    os << "Equilibrium did not converge: " << e.what() << endl;
    Messages::error(os.str());
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
  get_my_options().coupling = coupling;

  solveropts.set_option("nonlin_max_it", max_it);
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
          get_physical_model(elem->subdomain_id()));

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
DriftDiffusion::do_set_to_remembered_solution(ID id)
{
  // call the default implementation
  SimulationInterface::do_set_to_remembered_solution(id);

  if (do_local_scaling_)
    build_local_scaling();

}

void
DriftDiffusion::do_print_info(void)
{
  parse_const_options();
  parse_options();

  Options& myopts = get_my_options();

  Messages::newline();

  ostringstream os;
  os << "solving for : ";
  if (myopts.coupling & POISSON)
    os << "poisson ";
  if (myopts.coupling & ELECTRONS)
    os << "electrons ";
  if (myopts.coupling & HOLES)
    os << "holes ";

  os << endl;

  if (do_local_scaling_)
    os << "using local scaling";

  Messages::info(os.str());

}





void
DriftDiffusion::parse_const_options(void)
{

  const ModelOptions& opts = get_options();
  Options& myopts = get_my_options();

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

  string qrule = get_mesh().mesh_dimension() == 1 ? "trapez" : "gauss";
  qrule = opts.get_option("quadrature_rule", qrule);
  if (qrule == "gauss")
    myopts.quadrature_type = QGAUSS;
  else if (qrule == "trapez")
    myopts.quadrature_type = QTRAP;
  else
    throw InitFailedException("Unknown quadrature rule");

}



void
DriftDiffusion::parse_options(void)
{
  PerfLog perf_log("DriftDiffusion parse_options()", false);
  perf_log.start_event("parse");


  const ModelOptions& opts = get_options();
  Options& myopts = get_my_options();

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", 5));

  string coupling = opts.get_option("coupling", "");
  if (coupling == "full")
    myopts.coupling = FULLYCOUPLED;
  else if (coupling == "poisson")
    myopts.coupling = POISSON;
  else if (coupling == "electrons")
  {
    myopts.coupling = ECURRENT | POISSON;
    if (opts.get_option("local_equilibrium", false))
      _useparticle = 'e';
  }
  else if (coupling == "holes")
  {
    myopts.coupling = HCURRENT | POISSON;
    if (opts.get_option("local_equilibrium", false))
      _useparticle = 'h';
  }
  else if (coupling == "current")
    myopts.coupling = CURRENTS;


  myopts.mesh_refinement = opts.get_option("mesh_refinement",
      myopts.mesh_refinement);

  myopts.exact_newton = opts.get_option("exact_newton", myopts.exact_newton);

  perf_log.stop_event("parse");
}





void
DriftDiffusion::rebuild_equation_system(void)
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
  system.add_variable("fermi_e", libMeshEnums::FIRST);
  system.add_variable("fermi_h", libMeshEnums::FIRST);

  system.add_vector("old_sol");

  // finally initialize the newly created system
  system.init();


  _rebuild_eq_system = false;

}




void
DriftDiffusion::do_init(void)
{

  _device = &get_environment().get_device();

  parse_const_options();

  rebuild_equation_system();

  find_dirichlet_nodes();
  find_dielectric_boundary_nodes();

  set<const Boundary*> real_contacts;
  // prepare the _boundary_currents
  // we will rely on the fact that it contains an entry for every boundary
  // later on !!!
  {
    SimulationEnvironment::BoundaryIterator
      it(get_environment().boundaries_begin());
    const SimulationEnvironment::BoundaryIterator
      end(get_environment().boundaries_end());
    for ( ; it != end; ++it)
    {
      BoundaryProperties* bd = it->second->get_boundary_properties(get_id());
      if (bd != NULL)
      {
        ElectricalContact* contact = dynamic_cast<ElectricalContact*>(bd);
#ifndef DEBUG
        if (contact->is_real_contact())
#endif
        {
          _boundary_currents[it->second] = 0.0;
          _voltages[it->second] = 0.0;
          real_contacts.insert(it->second);
        }
      }
    }
  }

  get_environment().update_boundary_element_map(real_contacts);

}


void
DriftDiffusion::do_setup_solution_variables(void)
{
  // declare solution variables
  declare_solution(ElPotential, REAL, NODES, "V");
  declare_solution(eQFermi, REAL, NODES, "eV");
  declare_solution(hQFermi, REAL, NODES, "eV");
  declare_solution(ElField, VECTOR, CELL, "V/cm");
  declare_solution(Eg, REAL, NODES, "eV");
  declare_solution(Ec, REAL, NODES, "eV");
  declare_solution(Ev, REAL, NODES, "eV");
  declare_solution(Ec0, REAL, NODES, "eV");
  declare_solution(Ev0, REAL, NODES, "eV");

  declare_solution(Polarization, VECTOR, CELL, "C/m^2");

  declare_solution(eDensity, REAL, NODES, "cm^-3");
  declare_solution(hDensity, REAL, NODES, "cm^-3");

  declare_solution(eMobility, REAL, NODES, "cm^2/(V*s)");
  declare_solution(hMobility, REAL, NODES, "cm^2/(V*s)");

  declare_solution(eConductivity, REAL, NODES, "S/cm");
  declare_solution(hConductivity, REAL, NODES, "S/cm");

  declare_solution(CurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(eCurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(hCurrentDensity, VECTOR, CELL, "A/cm^2");
  //declare_solution(eFlux, VECTOR, CELL, "1/(s*cm^2)");
  //declare_solution(hFlux, VECTOR, CELL, "1/(s*cm^2)");

  declare_solution(IonizedDonors, REAL, NODES, "cm^-3");
  declare_solution(IonizedAcceptors, REAL, NODES, "cm^-3");

  declare_solution(eThElPower, REAL, NODES, "eV/K");
  add_alias("Pn", eThElPower);
  declare_solution(hThElPower, REAL, NODES, "eV/K");
  add_alias("Pp", hThElPower);

  declare_solution(eJoule, REAL, NODES, "W/cm^3");
  declare_solution(hJoule, REAL, NODES, "W/cm^3");

  declare_solution(ePowerFlux, VECTOR, NODES, "W/cm^2");
  declare_solution(hPowerFlux, VECTOR, NODES, "W/cm^2");

  declare_solution(ePeltier, REAL, NODES, "W/cm^3");
  declare_solution(hPeltier, REAL, NODES, "W/cm^3");

  declare_solution(RecombHeat, REAL, NODES, "W/cm^3");

  declare_solution(NetRecombination, REAL, NODES, "1/(s*cm^3)");
  // add the single recombination rates
  {
    bool plot_rec = plot_solution(NetRecombination);
    const set<PhysicalModel*>& pm = get_physical_models();
    set<PhysicalModel*>::const_iterator it(pm.begin());
    set<PhysicalModel*>::const_iterator end(pm.end());
    for ( ; it != end; ++it)
    {
      DriftDiffusionProperties* sc =
          static_cast<DriftDiffusionProperties*>(*it);

      vector<ID> ids;
      int n = sc->get_net_recombination_rate_IDs(ids);

      for (int i = 0; i < n; i++)
      {
        ID id = static_cast<ID>(NetRecombination) + ids[i];
        const std::string& name =
            sc->get_recombination_model(ids[i])->get_default_name();
        // if recombination should be plotted, add it also to the plot variables
        if (plot_rec) add_plot_variable(name);
        declare_solution_ext(name, id, SolutionDescriptor::REAL,
            SolutionDescriptor::NODES, "1/(s*cm^3)");
        _recombination_ids.insert(id);
      }
    }
  }

  // the contact currents/voltages
  {
    unsigned int dim = get_mesh().mesh_dimension();
    string units("A");
    switch (dim)
    {
      case 1:
        units = "A/cm^2";
        break;
      case 2:
        if (get_environment().get_device().get_symmetry()
            != TiberCad::CYLINDRICAL)
          units = "A/cm";
        break;
      default:
        break;
    }

    declare_solution(ContactCurrent, REAL, GLOBAL, units);
    add_alias("ContactCurrents", ContactCurrent);
    bool plot_curr = plot_solution(ContactCurrent);

    SimulationEnvironment::BoundaryIterator it(get_environment().boundaries_begin());
    const SimulationEnvironment::BoundaryIterator end(get_environment().boundaries_end());
    unsigned int i = 1;

    for ( ; it != end; ++it)
    {
      BoundaryProperties* bd = it->second->get_boundary_properties(get_id());
      if (bd != NULL)
      {
        ElectricalContact* contact = dynamic_cast<ElectricalContact*>(bd);
        if (contact->is_real_contact())
        {
          ID id = static_cast<ID>(ContactCurrent) + i;
          string name(it->second->get_name() + ".current");
          // if currents should be plotted, add it also to the plot variables
          if (plot_curr) add_plot_variable(name);
          declare_solution_ext(name, id, SolutionDescriptor::REAL,
              SolutionDescriptor::GLOBAL, units);
          i++;
        }
      }
    }
  }
}






void
DriftDiffusion::do_newton(void)
{

  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system =
    es.get_system<TiberNonlinearSystem>(get_equation_system_name());


  system.set_options(get_solver_options());
  system.solve();
}








void
DriftDiffusion::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& points)
{

  unsigned int np = points.size();

  TiberNonlinearSystem* system;
  system = &get_equation_systems().get_system<TiberNonlinearSystem>(
      get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  DriftDiffusionProperties* sc =
    dynamic_cast<DriftDiffusionProperties*>(
        get_physical_model(subdomain));

  assert(sc != NULL);

  sc->reinit(elem);

  fe->reinit(elem, &points);


  vector<double> T_nodes = sc->get_temperature_at_nodes();

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_ep, ep_var);

  const unsigned int n_dofs = dof_indices_u.size();

  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();
  double vol0 = get_scaling().get_calc_mesh_units();
  switch (dim)
  {
    case 3:
      vol0 *= get_scaling().get_calc_mesh_units();
    case 2:
      vol0 *= get_scaling().get_calc_mesh_units();
      break;
  }

  // cell data variables (to be integrated)
  RealGradient jn(0);
  RealGradient jp(0);
  RealGradient el_field(0);
  RealVectorValue polariz(0);

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    double en = 0.0;
    double ep = 0.0;
    double T  = 0.0;
    RealGradient e_field(0);
    RealGradient grad_en_loc(0);
    RealGradient grad_ep_loc(0);
    RealGradient grad_T_loc(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices_u[i]);
      en += phi[i][n] * solution(dof_indices_en[i]);
      ep += phi[i][n] * solution(dof_indices_ep[i]);

      grad_en_loc += dphi[i][n] * solution(dof_indices_en[i]);
      grad_ep_loc += dphi[i][n] * solution(dof_indices_ep[i]);

      e_field += dphi[i][n] * solution(dof_indices_u[i]);

      grad_T_loc += dphi[i][n] * T_nodes[i];

      T +=  phi[i][n] * T_nodes[i];

    }

    // scale the potential back
    u *= phi0;
    en *= phi0;
    ep *= phi0;
    e_field *= -phi0;
    grad_en_loc *= phi0;
    grad_ep_loc *= phi0;


    sc->set_coordinates(real_pts[n]);

    sc->set_potentials(u, en, ep);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_e(grad_en_loc);
    sc->set_grad_fermi_h(grad_ep_loc);

    sc->calculate_densities();

    double edens = (sc->is_dielectric() ? 0.0 : sc->get_electron_density());
    double hdens = (sc->is_dielectric() ? 0.0 : sc->get_hole_density());

    sc->calculate_mobilities();

    sc->compute_thermoelectric_powers();
    double Pn =  sc->get_electron_thermoelectric_power();
    double Pp =  sc->get_hole_thermoelectric_power();

    sc->calculate_net_recombination_rates();

    double sigma_e = Constants::e * edens * sc->get_electron_mobility();
    double sigma_h = Constants::e * hdens * sc->get_hole_mobility();

    RealGradient dfn = grad_en_loc + Pn * grad_T_loc;
    RealGradient dfp = grad_ep_loc + Pp * grad_T_loc;

    RealGradient jn_loc = -sigma_e * dfn;
    RealGradient jp_loc = -sigma_h * dfp;
    jn += jn_loc;
    jp += jp_loc;

    el_field += e_field;
    polariz += sc->get_total_polarization();


    if (values.count(ElPotential))
      values[ElPotential][n] = u;

    if (values.count(eQFermi))
      values[eQFermi][n] = -en;

    if (values.count(hQFermi))
      values[hQFermi][n] = -ep;

    if (values.count(Ec))
      values[Ec][n] = sc->get_conduction_band_edge() - u;

    if (values.count(Ev))
      values[Ev][n] = sc->get_valence_band_edge() - u;

    if (values.count(Ec0))
      values[Ec0][n] = sc->get_conduction_band_edge();

    if (values.count(Ev0))
      values[Ev0][n] = sc->get_valence_band_edge();

    if (values.count(Eg))
      values[Eg][n] =
        sc->get_conduction_band_edge() - sc->get_valence_band_edge();
/*
    if (ids.count(VBANDEDGEINTR))
      values[n][VBANDEDGEINTR] = sc->get_valence_band_edge();

    if (ids.count(CBANDEDGEINTR))
      values[n][CBANDEDGEINTR] = sc->get_conduction_band_edge();
*/
    if (values.count(eDensity))
      values[eDensity][n] = edens;

    if (values.count(hDensity))
      values[hDensity][n] = hdens;

    if (values.count(eMobility))
      values[eMobility][n] = sc->get_electron_mobility();

    if (values.count(hMobility))
      values[hMobility][n] = sc->get_hole_mobility();

    if (values.count(eConductivity))
      values[eConductivity][n] = sigma_e;

    if (values.count(hConductivity))
      values[hConductivity][n] = sigma_h;

    bool ionized_donors = values.count(IonizedDonors);
    bool ionized_acceptors = values.count(IonizedAcceptors);
    if (ionized_donors || ionized_acceptors)
    {
      sc->calculate_ionized_dopants();
      if (ionized_donors)
        values[IonizedDonors][n] = sc->get_ionized_donor_density();

      if (ionized_acceptors)
        values[IonizedAcceptors][n] = sc->get_ionized_acceptor_density();
    }

    if (values.count(eJoule))
      values[eJoule][n] = -(jn_loc * dfn);

    if (values.count(hJoule))
      values[hJoule][n] = -(jp_loc * dfp);

    if (values.count(ePowerFlux))
    {
      values[ePowerFlux][3 * n] = (Pn * T + en) * jn_loc(0);
      values[ePowerFlux][3 * n + 1] = (Pn * T + en) * jn_loc(1);
      values[ePowerFlux][3 * n + 2] = (Pn * T + en) * jn_loc(2);
    }

    if (values.count(hPowerFlux))
    {
      values[hPowerFlux][3 * n] = (Pp * T + ep) * jp_loc(0);
      values[hPowerFlux][3 * n + 1] = (Pp * T + ep) * jp_loc(1);
      values[hPowerFlux][3 * n + 2] = (Pp * T + ep) * jp_loc(2);
    }

    if (values.count(eThElPower))
      values[eThElPower][n] = Pn;

    if (values.count(hThElPower))
      values[hThElPower][n] = Pp;

    {
      bool get_recomb = values.count(NetRecombination);
      double tot_rec = 0;

      // loop over all recombination ids
      set<ID>::const_iterator rec_it(_recombination_ids.begin());
      for ( ; rec_it != _recombination_ids.end(); ++rec_it)
      {
        bool requested = values.count(*rec_it);
        double rec = 0;
        if (get_recomb || requested)
          rec = sc->get_net_recombination_rate(*rec_it - NetRecombination);

        if (requested)
          values[*rec_it][n] = rec;

        if (get_recomb)
          tot_rec += rec;
      }

      if (get_recomb)
        values[NetRecombination][n] = tot_rec;
    }


    {
      bool ept = values.count(ePeltier);
      bool hpt = values.count(hPeltier);

      if (ept || hpt)
        sc->compute_thermoelectric_power_gradient();

      if (ept)
      {
        RealGradient PnGrad = sc->get_electron_thermoelectric_power_gradient();
        values[ePeltier][n] =  -T * (PnGrad * jn_loc);
      }

      if (hpt)
      {
        RealGradient PpGrad = sc->get_hole_thermoelectric_power_gradient();
        values[hPeltier][n] =  -T * (PpGrad * jp_loc);
      }
    }

    if (values.count(RecombHeat))
    {
      vector<ID> rec_model_ids;
      int n_rec = sc->get_net_recombination_rate_IDs(rec_model_ids);
      double rec = 0.0;
      for (int i = 0; i < n_rec; i++)
	rec += sc->get_net_recombination_rate(rec_model_ids[i]);

      values[RecombHeat][n] = Constants::e * rec * (ep - en + T * (Pp - Pn));
    }

  }


  // the cell based solutions
  // for now, we make a mean value

  if (values.count(ElField))
  {
    values[ElField][0] = el_field(0) / np;
    values[ElField][1] = el_field(1) / np;
    values[ElField][2] = el_field(2) / np;
  }

  if (values.count(Polarization))
  {
    values[Polarization][0] = polariz(0) / np;
    values[Polarization][1] = polariz(1) / np;
    values[Polarization][2] = polariz(2) / np;
  }

  if (values.count(CurrentDensity))
  {
    values[CurrentDensity][0] = (jn(0) + jp(0)) / np;
    values[CurrentDensity][1] = (jn(1) + jp(1)) / np;
    values[CurrentDensity][2] = (jn(2) + jp(2)) / np;
  }

  if (values.count(eCurrentDensity))
  {
    values[eCurrentDensity][0] = jn(0) / np;
    values[eCurrentDensity][1] = jn(1) / np;
    values[eCurrentDensity][2] = jn(2) / np;
  }

  if (values.count(hCurrentDensity))
  {
    values[hCurrentDensity][0] = jp(0) / np;
    values[hCurrentDensity][1] = jp(1) / np;
    values[hCurrentDensity][2] = jp(2) / np;
  }

}







void
DriftDiffusion::calculate_currents_rstf(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;

  // reset currents
  {
    ContactData::iterator it =
      _boundary_currents.begin();
    for ( ; it != _boundary_currents.end(); ++it)
      (*it).second = 0.0;
  }

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  SimulationEnvironment& env = get_environment();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();


  const double phi0 = get_scaling().get_potential_scaling();


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  FEType fe_type = system->variable_type(u_var);

  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
  AutoPtr<QBase> qrule(QBase::build(
        get_my_options().quadrature_type, dim, get_my_options().integration_order));
  fe->attach_quadrature_rule(qrule.get());


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


  BoundaryElementMap::iterator el(env.boundary_elements_begin());
  BoundaryElementMap::iterator end_el(env.boundary_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = elem->top_parent();

    const Boundary* boundary = el.get_boundary();

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));

    assert(sc != NULL);


    // in a dielectric we have no current
    if (sc->is_dielectric())
      continue;


    fe->reinit(elem);

    sc->reinit(elem);

    //Get the temperature given the element
    vector<double> T_nodes =  sc->get_temperature_at_nodes();

    // find the weight for each element node
    vector<double> weight(elem->n_nodes(), 0.0);
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
      if (env.is_node_on_boundary(elem->get_node(n), boundary))
        weight[n] = 1.0;


    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
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

      //Get the thermoelectric power
      sc->compute_thermoelectric_powers();
      double Pn =  sc->get_electron_thermoelectric_power() / phi0;
      double Pp =  sc->get_hole_thermoelectric_power() / phi0;

      // we put the minus here for convenience
      double sigma_e = -Constants::e * sc->get_electron_density() *
        sc->get_electron_mobility();
      double sigma_h = -Constants::e * sc->get_hole_density() *
        sc->get_hole_mobility();

      RealGradient je(JxW[qp] * phi0 * (sigma_e * (dEfn + Pn * dT)));
      RealGradient jh(JxW[qp] * phi0 * (sigma_h * (dEfp + Pp * dT)));

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        _boundary_currents[boundary] += (je + jh) * dphi[n][qp] * weight[n];

    } // end loop over quadrature points
  } // end loop over elements

}



void
DriftDiffusion::calculate_field_emission(void)
{

  ContactData fe_currents;
  ContactData::iterator it = _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
  {
    ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
        (*it).first->get_boundary_properties(get_id()));
    if (contact->has_field_emission())
      fe_currents[(*it).first] = 0.0;
  }

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = get_scaling().get_potential_scaling();


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");

  FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = get_my_options().integration_order;

  AutoPtr<QBase> qface(QBase::build(
        get_my_options().quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


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

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));

    assert(sc != NULL);

    /* we calculate field emission only in dielectrics */
    if (!sc->is_dielectric())
      continue;

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(top_parent, s);

      if (env.is_boundary(side))
      {

        Boundary* boundary = env.get_boundary(side);
        if (boundary == NULL)
          continue;

        ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
            boundary->get_boundary_properties(get_id()));

        // check if we should do something
        if (!contact->has_field_emission())
          break;

        sc->reinit(elem);

        fe_face->reinit(elem, s);

        int phi_size = phi.size();

        double current = 0.0;

        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {
          // get the solution value at the quadrature point
          double e_field = 0.0;
          for (unsigned int i = 0; i < phi_size; i++)
          {
            double tmp = dphi[i][qp] * face_normals[qp];
            e_field += tmp * solution(dof_indices_u[i]);
          }

          double F = phi0 * e_field;
          current += JxW[qp] * contact->calculate_field_emission(F);
        } // end loop over quadrature points

        fe_currents[boundary] += current;
      }
    } // end loop over elem sides
  } // end loop over elements


  for (it = fe_currents.begin(); it != fe_currents.end(); ++it)
  {
    ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
        (*it).first->get_boundary_properties(get_id()));

    contact->set_field_emission_current((*it).second);
  }
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
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = get_scaling().get_potential_scaling();


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = get_my_options().integration_order;

  AutoPtr<QBase> qface(QBase::build(
        get_my_options().quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


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


  BoundaryElementMap::iterator el(env.boundary_elements_begin());
  BoundaryElementMap::iterator end_el(env.boundary_elements_end());

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
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));

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

        //Get the temperature given the element
        vector<double> T_nodes = sc->get_temperature_at_nodes();

        // only for dim > 1 we need to integrate
        if (dim > 1)
        {
          fe_face->reinit(elem, s);

          int phi_size = phi.size();

          double current = 0.0;

          for (unsigned int qp = 0; qp < qface->n_points(); qp++)
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
            sc->set_grad_fermi_e(phi0 * dEfn);
            sc->set_grad_fermi_h(phi0 * dEfp);

            sc->calculate_densities();

            sc->calculate_mobilities();

	    //Get the thermoelectric power------------
	    sc->compute_thermoelectric_powers();
	    double Pn = sc->get_electron_thermoelectric_power() / phi0;
	    double Pp = sc->get_hole_thermoelectric_power() / phi0;

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
          sc->set_grad_fermi_e(phi0 * RealGradient(dEfn, 0, 0));
          sc->set_grad_fermi_h(phi0 * RealGradient(dEfp, 0, 0));

          sc->calculate_densities();

          sc->calculate_mobilities();

	  //Get the thermoelectric power------------
	  sc->compute_thermoelectric_powers();
	  double Pn = sc->get_electron_thermoelectric_power() / phi0;
	  double Pp = sc->get_hole_thermoelectric_power() / phi0;


          Real cond_e = Constants::e * sc->get_electron_mobility() *
            sc->get_electron_density();
          Real cond_h = Constants::e * sc->get_hole_mobility() *
            sc->get_hole_density();

          if (boundary != NULL)
          {
            //ElectricalContact* contact = dynamic_cast<ElectricalContact*>(
            //   boundary->get_boundary_properties(get_id()));
            _boundary_currents[boundary] = -phi0 *
              (cond_e * (dEfn + Pn * dT) + cond_h * (dEfp + Pp * dT));
          }
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
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const Options& params = get_my_options();

  // the scaling parameters to scale back the result
  const Scaling& scaling = get_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double x0 = scaling.get_length_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;


  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  AutoPtr<QBase> qrule(QBase::build(
        params.quadrature_type, dim, params.integration_order));
  fe->attach_quadrature_rule(qrule.get());


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
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));
    assert(sc != NULL);

    sc->reinit(elem);

    fe->reinit(elem);

    assert(elem->n_nodes() == dof_indices_u.size());

    RealGradient field(0.0);
    for (unsigned int i = 0; i < dof_indices_u.size(); i++)
      field += dphi[i][0] * solution(dof_indices_u[i]);
    field *= -phi0;


    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
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


      const double penalty_value = 1e56;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        if (_dirichlet_dofs.count(dof_indices_en[i]))
          local_scaling_[elem->get_node(i)][0] += penalty_value;
        else
          local_scaling_[elem->get_node(i)][0] +=
            sigma_e * (dphi[i][qp] * dphi[i][qp]);

        if (_dirichlet_dofs.count(dof_indices_ep[i]))
          local_scaling_[elem->get_node(i)][1] += penalty_value;
        else
        local_scaling_[elem->get_node(i)][1] +=
          sigma_h * (dphi[i][qp] * dphi[i][qp]);

        if (_dirichlet_dofs.count(dof_indices_u[i]))
          local_scaling_[elem->get_node(i)][2] += penalty_value;
        else
        local_scaling_[elem->get_node(i)][2] +=
          l2_eps * (dphi[i][qp] * dphi[i][qp]) -
          drho * phi[i][qp] * phi[i][qp];
      }


    } // end loop over quadrature points
  } // end loop over elements
}



NumericVector<double>&
DriftDiffusion::do_get_solution_vector(void)
{
  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  system->get_solution_vector().close();
  return system->get_solution_vector();
}


/*
// implementation taken from libmesh equation_systems.C
void
DriftDiffusion::build_nodal_results(const set<string>& variables,
    vector<double>& results, vector<string>& legend)
{

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = get_solution_vector();
  const NumericVector<Number>& oldsolution = system->get_vector("old_sol");

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  // TODO if some elements were coarsened, does this still work??
  const unsigned int nn  = mesh.n_nodes();

  legend.resize(variables.size());

  bool found_variable = false;

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

    found_variable = true;
  }

  int Ev = -1;
  if (variables.find("Ev") != varend)
  {
    Ev = n_vars;
    legend[n_vars] = "Ev";
    n_vars++;

    found_variable = true;
  }

  int Ec0 = -1;
  if (variables.find("Ec0") != varend)
  {
    Ec0 = n_vars;
    legend[n_vars] = "Ec0";
    n_vars++;

    found_variable = true;
  }

  int Ev0 = -1;
  if (variables.find("Ev0") != varend)
  {
    Ev0 = n_vars;
    legend[n_vars] = "Ev0";
    n_vars++;

    found_variable = true;
  }


  int Efn = -1;
  if (variables.find("QFermi_e") != varend)
  {
    Efn = n_vars;
    legend[n_vars] = "QFermi_e";
    n_vars++;

    found_variable = true;
  }

  int Efp = -1;
  if (variables.find("QFermi_h") != varend)
  {
    Efp = n_vars;
    legend[n_vars] = "QFermi_h";
    n_vars++;

    found_variable = true;
  }

  int phi = -1;
  if (variables.find("ElPotential") != varend)
  {
    phi = n_vars;
    legend[n_vars] = "electric_potential";
    n_vars++;

    found_variable = true;
  }

  int Eg = -1;
  if (variables.find("Eg") != varend)
  {
    Eg = n_vars;
    legend[n_vars] = "Eg";
    n_vars++;

    found_variable = true;
  }

  int edens = -1;
  if (variables.find("eDensity") != varend)
  {
    edens = n_vars;
    legend[n_vars] = "electron_density";
    n_vars++;

    found_variable = true;
  }

  int hdens = -1;
  if (variables.find("hDensity") != varend)
  {
    hdens = n_vars;
    legend[n_vars] = "hole_density";
    n_vars++;

    found_variable = true;
  }

  int Nd = -1;
  if (variables.find("Nd") != varend)
  {
    Nd = n_vars;
    legend[n_vars] = "ionized_donors";
    n_vars++;

    found_variable = true;
  }

  int Na = -1;
  if (variables.find("Na") != varend)
  {
    Na = n_vars;
    legend[n_vars] = "ionized_acceptors";
    n_vars++;

    found_variable = true;
  }

  int Pn = -1;
  if (variables.find("Pn") != varend)
  {
    Pn = n_vars;
    legend[n_vars] = "Pn";
    n_vars++;

    found_variable = true;
  }

  int Pp = -1;
  if (variables.find("Pp") != varend)
  {
    Pp = n_vars;
    legend[n_vars] = "Pp";
    n_vars++;

    found_variable = true;
  }

  int rho = -1;
  if (variables.find("charge_density") != varend)
  {
    rho = n_vars;
    legend[n_vars] = "total_charge_densitity";
    n_vars++;

    found_variable = true;
  }

  int rec = -1;
  int num_rec = 0;
  map<ID, string> rec_model_ids;
  if (variables.find("NetRecombination") != varend)
  {

    // look for all recombination models
    const set<PhysicalModel*>& pm = get_physical_models();

    set<PhysicalModel*>::const_iterator it(pm.begin());
    set<PhysicalModel*>::const_iterator end(pm.end());
    for ( ; it != end; ++it)
    {
      DriftDiffusionProperties* sc =
        static_cast<DriftDiffusionProperties*>(*it);

        vector<ID> ids;

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

    found_variable = true;
  }


  int mun = -1;
  if (variables.find("eMob") != varend)
  {
    mun = n_vars;
    legend[n_vars] = "electron_mobility";
    n_vars++;

    found_variable = true;
  }

  int mup = -1;
  if (variables.find("hMob") != varend)
  {
    mup = n_vars;
    legend[n_vars] = "hole_mobility";
    n_vars++;

    found_variable = true;
  }


  legend.resize(n_vars);

  results.resize(nn * n_vars);

  // We return immediately if there is nothing to plot
  if (!found_variable)
    return;

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
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

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
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));
      assert(sc != NULL);

      sc->reinit(elem);

      fe->reinit(elem);

      assert(elem->n_nodes() == dof_indices_u.size());

      RealGradient field(0.0);
      RealGradient grad_Fe(0.0);
      RealGradient grad_Fh(0.0);
      for (unsigned int i = 0; i < dof_indices_u.size(); i++)
      {
        field += dphi[i][0] * solution(dof_indices_u[i]);
        grad_Fe += dphi[i][0] * solution(dof_indices_en[i]);
        grad_Fh += dphi[i][0] * solution(dof_indices_ep[i]);
      }
      field *= -phi0;
      grad_Fe *= -phi0;
      grad_Fh *= -phi0;


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        double u  = phi0 * solution(dof_indices_u[n]);
        double en = phi0 * solution(dof_indices_en[n]);
        double ep = phi0 * solution(dof_indices_ep[n]);
        double oldu  = phi0 * oldsolution(dof_indices_u[n]);
        double olden = phi0 * oldsolution(dof_indices_en[n]);
        double oldep = phi0 * oldsolution(dof_indices_ep[n]);


        // prepare for calculating local properties
        sc->set_coordinates(elem->point(n));


        sc->set_potentials(u, en, ep);
        sc->set_old_potentials(oldu, olden, oldep);
        sc->set_electric_field(field);
        sc->set_grad_fermi_e(grad_Fe);
        sc->set_grad_fermi_h(grad_Fh);

        sc->calculate_densities();
        sc->compute_thermoelectric_powers();
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

	if (Pn != -1)
        {
          double nodal_val = sc->get_electron_thermoelectric_power();
          local[id + Pn] += nodal_val / conn;
        }

	if (Pp != -1)
        {
          double nodal_val = sc->get_hole_thermoelectric_power();
          local[id + Pp] += nodal_val / conn;
        }

        if (phi != -1)
        {
          local[id + phi] += u / conn;
        }

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
*/





/*
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
  const NumericVector<Number>& oldsolution = system->get_vector("old_sol");

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_active_elem();

  legend.resize(variables.size());

  // for each possible variable we set the vector index
  // -1 means, the variable should not be plotted
  unsigned int n_vars = 0;
  const set<string>::const_iterator varend(variables.end());

  bool found_variable = false;

  int BandEdges = -1;
  if (variables.find("BandEdges") != varend)
  {
    const Elem* elem = *mesh.active_elements_begin();

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(
          get_physical_model(elem->subdomain_id()));

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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
  }


  int J = -1;
  if ((variables.find("CurrentDensity") != varend) ||
      (variables.find("Current") != varend))
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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
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

    found_variable = true;
  }


  int PDens = -1;
  //if (variables.find("PowerDensity") != varend)
  //{
  //  PDens = n_vars;
  //  legend[n_vars] = "power_density[W*cm^-3]";
  //  n_vars++;
  //}


  legend.resize(n_vars);

  results.resize(nn * n_vars);

  // We return immediately if there are no variables to be plotted
  if (!found_variable)
    return;


  // the scaling parameters to scale back the result
  double phi0 = get_scaling().get_potential_scaling();

  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

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
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));

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
      oldu  += phi[i][0] * oldsolution(dof_indices_u[i]);
      olden += phi[i][0] * oldsolution(dof_indices_en[i]);
      oldep += phi[i][0] * oldsolution(dof_indices_ep[i]);
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
}
*/


void
DriftDiffusion::build_integrated_quantities(vector<double>& values)
{

  if (plot_solution("ContactCurrent") || plot_solution("ContactCurrents")
      || plot_solution("current"))
  {
    calculate_currents();

    values.resize(_boundary_currents.size());

    // we need alphabetic order !!
    map<string, double> currs;

    ContactData::iterator it(_boundary_currents.begin());
    const ContactData::iterator end(_boundary_currents.end());
    for (; it != end; ++it)
      currs[it->first->get_name()] = it->second * it->first->get_area_factor();


    map<string, double>::iterator nameit(currs.begin());
    map<string, double>::iterator nameend(currs.end());
    for (unsigned int id = 0; nameit != nameend; ++nameit, id++)
      values[id] = nameit->second;

  }
}




void
DriftDiffusion::get_solution_secure(map<ID, vector<double> >& values)
{
  vector<string> tokens;

  map<ID, vector<double> >::iterator mapit(values.begin());
  const map<ID, vector<double> >::iterator mapend(values.end());
  for ( ; mapit != mapend; ++mapit)
  {
    ID id = mapit->first;
    const SolutionDescriptor& descr = get_solution_descriptor(id);
    Utils::tokenize(descr.name(), tokens);

    ContactData::iterator it(_boundary_currents.begin());
    const ContactData::iterator end(_boundary_currents.end());
    for (; it != end; ++it)
    {
      if (tokens[0] == it->first->get_name())
      {
        double curr = it->second * it->first->get_area_factor();
        values[id] = vector<double>(1, curr);
        break;
      }
    }
  }
}



void
DriftDiffusion::calculate_currents(void)
{
  if (get_my_options().current_calculation == RSTF)
    calculate_currents_rstf();
  else
    calculate_currents_surfint();

  calculate_field_emission();
}


void
DriftDiffusion::build_integrated_quantities_description(
    vector<string>& legend,
    vector<string>& description)
{

  if (plot_solution("ContactCurrent") || plot_solution("ContactCurrents")
      || plot_solution("current"))
  {
    legend.resize(_boundary_currents.size());

    // we make first a set to order the contacts alphabetically
    set<string> bds;

    ContactData::iterator it(_boundary_currents.begin());
    const ContactData::iterator end(_boundary_currents.end());
    for (; it != end; ++it)
      bds.insert(it->first->get_name());

    set<string>::iterator nameit(bds.begin());
    set<string>::iterator nameend(bds.end());
    for (unsigned int id = 0; nameit != nameend; ++nameit, id++)
      legend[id] = *nameit;

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
        if (get_environment().get_device().get_symmetry()
            != TiberCad::CYLINDRICAL)
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
  const MeshBase& mesh = get_mesh();

  TiberNonlinearSystem& system = es.get_system<TiberNonlinearSystem>(
      get_equation_system_name());


  // the current solutions
  NumericVector<Number>& solution = system.get_solution_vector();

  DofList& dirichlet_dofs = _dirichlet_dofs;

  DofList::const_iterator dof_it;
  const DofList::const_iterator dof_end =
    dirichlet_dofs.end();


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
        bool is_dirichlet_node =
            (dirichlet_dofs.find(dof_indices_u[i]) != dof_end) ? true : false;
        is_dirichlet_node |= (dirichlet_dofs.find(dof_indices_n[i]) != dof_end) ? true : false;
        is_dirichlet_node |= (dirichlet_dofs.find(dof_indices_p[i]) != dof_end) ? true : false;
        if (is_dirichlet_node)
        {
          Boundary* bd = get_environment().get_boundary(elem->get_node(i));
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
  const MeshBase& mesh = get_mesh();
  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = static_cast<TiberNonlinearSystem&>(
      eq_sys.get_system(get_equation_system_name()));

  const unsigned int dim = mesh.mesh_dimension();

  const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  const Options& params = get_my_options();
  Options& options = get_my_options();


  const NumericVector<Number>& oldx = system.get_vector("old_sol");

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
  // x 1e4 because we calculate in cm, but P comes in C/m^2
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
  AutoPtr<QBase> qrule(QBase::build(
        params.quadrature_type, dim, integration_order));
  fe->attach_quadrature_rule(qrule.get());

  // the finite element for boundary integration
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));

  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;

  AutoPtr<QBase> qface(QBase::build(
        params.quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


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
  DenseVector<Number> oldX;

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
    oldXu(oldX),
    oldXn(oldX),
    oldXp(oldX);


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
    oldX.resize(n_dofs_tot);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);
    dof_map.extract_local_vector(oldx, dof_indices, oldX);

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
    if (_useparticle == 'h')
      Xn.reposition(2 * n_dofs, n_dofs);
    else
      Xn.reposition(n_dofs, n_dofs);
    if (_useparticle == 'e')
      Xp.reposition(n_dofs, n_dofs);
    else
      Xp.reposition(2 * n_dofs, n_dofs);
    //
    oldXu.reposition(0, n_dofs);
    oldXn.reposition(n_dofs, n_dofs);
    oldXp.reposition(2 * n_dofs, n_dofs);



    DriftDiffusionProperties* sc =
      dynamic_cast<DriftDiffusionProperties*>(get_physical_model(subdomain));

    assert(sc != NULL);
    sc->reinit(elem);


    // Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_at_nodes();


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
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      Real oldu  = 0.0;
      Real olden = 0.0;
      Real oldep = 0.0;
      RealGradient e_field(0);
      RealGradient grad_en(0);
      RealGradient grad_ep(0);
      //RealGradient olde_field(0);
      //RealGradient oldgrad_en(0);
      //RealGradient oldgrad_ep(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * Xu(i);
        en += phi[i][qp] * Xn(i);
        ep += phi[i][qp] * Xp(i);
        oldu  += phi[i][qp] * oldXu(i);
        olden += phi[i][qp] * oldXn(i);
        oldep += phi[i][qp] * oldXp(i);
        e_field -= dphi[i][qp] * Xu(i);
        grad_en += dphi[i][qp] * Xn(i);
        grad_ep += dphi[i][qp] * Xp(i);
        //olde_field -= dphi[i][qp] * oldXu(i);
        //oldgrad_en += dphi[i][qp] * oldXn(i);
        //oldgrad_ep += dphi[i][qp] * oldXp(i);
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
      sc->set_old_potentials(phi0 * oldu, phi0 * olden, phi0 * oldep);

      double grad_fac = phi0 / x0;
      sc->set_electric_field(grad_fac * e_field);
      sc->set_grad_fermi_e(grad_fac * grad_en);
      sc->set_grad_fermi_h(grad_fac * grad_ep);

      sc->calculate_densities();

      long double n = sc->get_electron_density();
      long double p = sc->get_hole_density();
      //double Nd = sc->get_ionized_donor_density();
      //double Na = sc->get_ionized_acceptor_density();

      // calculate all local properties
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      // Get the thermoelectric power
      sc->compute_thermoelectric_powers();
      double eTEpower =  sc->get_electron_thermoelectric_power() / phi0;
      double hTEpower =  sc->get_hole_thermoelectric_power() / phi0;

      double epsilon = sc->get_relative_permittivity();
      long double l2_eps = l2 * epsilon;

      long double Rn = sc->get_net_electron_recombination_rate();
      //Rn = (fabs(Rn) < 1.0e-3) ? 0.0 : Rn;
      long double Rp = sc->get_net_hole_recombination_rate();
      //Rp = (fabs(Rp) < 1.0e-3) ? 0.0 : Rp;
      if (_useparticle != 'b')
        Rn = Rp = 0.0;


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

        if (_useparticle == 'e')
          drho[1] = -drho[0];
        else if (_useparticle == 'h')
          drho[2] = -drho[0];


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
        // the factor phi_0 comes from the derivative with respect to the potential
        long double dsigma_e = J * phi0 / (mu0 * C0_e) * mue * dn_dphi;
        long double dsigma_h = J * phi0 / (mu0 * C0_h) * muh * dp_dphi;

        // field dependent mobility
        // the factor phi_0 / x0 comes from the derivative with respect to the
        // gradient of the potential
        RealGradient dmu_e(0);
        RealGradient dmu_h(0);
        if (dim > 1)
        {
          sc->get_electron_mobility_derivatives(dmu_e);
          dmu_e *= J * phi0 / (mu0 * C0_e) * n / x0;
          sc->get_hole_mobility_derivatives(dmu_h);
          dmu_h *= J * phi0 / (mu0 * C0_h) * p / x0;
        }


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double lap_e = (dphi[i][qp] * grad_en) / local_scaling[i][0];
          double lap_h = (dphi[i][qp] * grad_ep) / local_scaling[i][1];
          long double dsigma_e_x_lap = dsigma_e * lap_e;
          long double dsigma_h_x_lap = dsigma_h * lap_h;

          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            //

            if (_options.exact_newton)
            {
              long double dsigma_e_x_phi = dsigma_e_x_lap * phi[j][qp];
              long double dsigma_h_x_phi = dsigma_h_x_lap * phi[j][qp];

              long double dmu_e_x_dphi = dmu_e * dphi[j][qp];
              long double dmu_h_x_dphi = dmu_h * dphi[j][qp];

              if (coupling & ECURRENT)
              {
                if (coupling & POISSON)
                  Knu(i,j) += dsigma_e_x_phi;

                Knn(i,j) += dmu_e_x_dphi * lap_e - dsigma_e_x_phi;
              }

              if (coupling & HCURRENT)
              {
                if (coupling & POISSON)
                  Kpu(i,j) += dsigma_h_x_phi;

                Kpp(i,j) += dmu_h_x_dphi * lap_h - dsigma_h_x_phi;
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
              if (coupling & POISSON)
                Knu(i,j) -= dRn[0] * phi_i_x_phi_j / local_scaling[i][0];

              Knn(i,j) -= dRn[1] * phi_i_x_phi_j / local_scaling[i][0];

              if (coupling & HCURRENT)
                Knp(i,j) -= dRn[2] * phi_i_x_phi_j / local_scaling[i][0];
            }

            if (coupling & HCURRENT)
            {
              if (coupling & POISSON)
                Kpu(i,j) += dRp[0] * phi_i_x_phi_j / local_scaling[i][1];

              if (coupling & ECURRENT)
                Kpn(i,j) += dRp[1] * phi_i_x_phi_j / local_scaling[i][1];

              Kpp(i,j) += dRp[2] * phi_i_x_phi_j / local_scaling[i][1];
            }

          }
        }

      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        long double J_x_rho;
        J_x_rho = J * sc->get_charge_density() / C0;

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
            if (_useparticle == 'h')
              Fn(i) -= 0;
            else
              Fn(i) -= Xn(i);

          if (coupling & HCURRENT)
            Fp(i) += net_recomb_h;
          else
            if (_useparticle == 'e')
              Fp(i) -= 0;
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

      // for the new implementation
      //PhysicalModel* sm = get_surface_model(elem, s);

      // is this a boundary?
      if (environment.is_boundary(side))
      {
        Boundary* boundary = environment.get_boundary(side);

        ElectricalContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<ElectricalContact*>(
              boundary->get_boundary_properties(get_id()));

        // we need to know if it is an outer boundary
        bool true_boundary = environment.is_outer_boundary(side);
        if (contact != NULL)
          contact->is_outer_boundary(true_boundary);

        //PhysicalModel* pm = get_surface_model(elem, s);

        fe_face->reinit(elem, s);

        // calculate the fluxes on the nodes
        if ((contact != NULL) && (dim > 1))
        {
          AutoPtr<Elem> side(elem->build_side(s));

          vector<Point> p(side->n_nodes());
          for (unsigned int i = 0; i < side->n_nodes(); i++)
            p[i] = side->point(i);

/*
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

            sc->compute_thermoelectric_powers();
            double Pn =  sc->get_electron_thermoelectric_power();
            double Pp =  sc->get_hole_thermoelectric_power();

            if (coupling & ECURRENT)
              nodal_flux_n[side->node(i)] =
                (sigma_e * grad_en + Pn * grad_T) * face_normals[0] / x0;
            if (coupling & HCURRENT)
              nodal_flux_p[side->node(i)] =
                -(sigma_h * grad_ep + Pp * grad_T) * face_normals[0] / x0;
          }
*/
        }


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
          for (unsigned int qp = 0; qp < qface->n_points(); qp++)
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
                e_field -= dphi_face[i][qp] * Xu(i);
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
              double Pn =  sc->get_electron_thermoelectric_power() / phi0;
              double Pp =  sc->get_hole_thermoelectric_power() / phi0;

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
                value[0] = c / (x0 * C0);
              }
              if (coupling & ECURRENT)
              {
                contact->get_normal_derivative(FERMIE, a, c);
                contact->get_derivatives_of_normal_derivative(FERMIE,
                    dcoeff[1], dvalue[1]);
                coeff[1] = a; // ???
                value[1] = c / (x0 * R0_e);
              }
              if (coupling & HCURRENT)
              {
                contact->get_normal_derivative(FERMIH, a, c);
                contact->get_derivatives_of_normal_derivative(FERMIH,
                    dcoeff[2], dvalue[2]);
                coeff[2] = a; // ???
                value[2] = c / (x0 * R0_h);
              }
            }



            // the jacobian x weight x scaling
            double J = JxW_face[qp];
/*
            // first the contributions to Ke_ij
            for (unsigned int i = 0; i < n_dofs; i++)
            {
              for (unsigned int j = 0; j < n_dofs; j++)
              {

                Real phi_i_x_phi_j =
                  J * phi_face[i][qp] * phi_face[j][qp];

                if (coupling & POISSON)
                  Kuu(i,j) += l2_eps * coeff[0] * phi_i_x_phi_j
                    / local_scaling[i][2];

                //if (coupling & ECURRENT)
                //  Knn(i,j) += coeff[1] * phi_i_x_phi_j / local_scaling[i][0];

                //if (coupling & HCURRENT)
                //  Kpp(i,j) += coeff[2] * phi_i_x_phi_j / local_scaling[i][1];
              }
            }
*/

            // contribution to the jacobian
            if (jacobian != NULL)
            {
              double scale_u = J * phi0 / x0 / C0;
              double scale_n = J * phi0 / (x0 * R0_e);
              double scale_p = J * phi0 / (x0 * R0_h);

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                double fac_u = scale_u / local_scaling[i][2];
                double fac_n = scale_n / local_scaling[i][0];
                double fac_p = scale_p / local_scaling[i][1];

                for (unsigned int j = 0; j < n_dofs; j++)
                {

                  Real phi_i_x_phi_j =
                    phi_face[i][qp] * phi_face[j][qp];

                  if (coupling & POISSON)
                  {
                    Kuu(i,j) -= fac_u * dvalue[0][0] * phi_i_x_phi_j;

                    if (coupling & ECURRENT)
                      Kun(i,j) -= fac_u * dvalue[0][1] * phi_i_x_phi_j;

                    if (coupling & HCURRENT)
                      Kup(i,j) -= fac_u * dvalue[0][2] * phi_i_x_phi_j;
                  }


                  ////
                  // NOTE:
                  //   the signs are inverted here because outflow means recombination
                  //   and the normal point outwards, giving a positiv current for
                  //   outflow
                  if (coupling & ECURRENT)
                  {
                    if (coupling & POISSON)
                      Knu(i,j) -= fac_n * dvalue[1][0] * phi_i_x_phi_j;

                    Knn(i,j) -= fac_n * dvalue[1][1] * phi_i_x_phi_j;

                    if (coupling & HCURRENT)
                      Knp(i,j) -= fac_n * dvalue[1][2] * phi_i_x_phi_j;
                  }


                  if (coupling & HCURRENT)
                  {
                    if (coupling & POISSON)
                      Kpu(i,j) += fac_p * dvalue[2][0] * phi_i_x_phi_j;

                    if (coupling & ECURRENT)
                      Kpn(i,j) += fac_p * dvalue[2][1] * phi_i_x_phi_j;

                    Kpp(i,j) += fac_p * dvalue[2][2] * phi_i_x_phi_j;
                  }

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
              double value_n = J * value[1];
              double value_p = J * value[2];

              for (unsigned int i = 0; i < n_dofs; i++)
              {
                if (coupling & POISSON)
                  Fu(i) -= value_u * phi_face[i][qp] / local_scaling[i][2];

                if (coupling & ECURRENT)
                  Fn(i) -= value_n * phi_face[i][qp] / local_scaling[i][0];

                if (coupling & HCURRENT)
                  Fp(i) += value_p * phi_face[i][qp] / local_scaling[i][1];
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
          sc->set_coordinates(elem->point(s));
          sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);

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
          sc->set_grad_fermi_e(phi0 / x0 * RealGradient(grad_en, 0.0, 0.0));
          sc->set_grad_fermi_h(phi0 / x0 * RealGradient(grad_ep, 0.0, 0.0));
          sc->calculate_densities();
          sc->calculate_mobilities();

          // we put the phi0 here for convenience
          double sigma_e = phi0 * sc->get_electron_mobility() *
            sc->get_electron_density();
          double sigma_h = phi0 * sc->get_hole_mobility() *
            sc->get_hole_density();

          sc->compute_thermoelectric_powers();
          double Pn =  sc->get_electron_thermoelectric_power() / phi0;
          double Pp =  sc->get_hole_thermoelectric_power() / phi0;

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

        if (contact != NULL)
        {
          if ((contact->get_type(POTENTIAL) == ElectricalContact::DIRICHLET)
              || (contact->get_type(FERMIE) == ElectricalContact::DIRICHLET)
              || (contact->get_type(FERMIH) == ElectricalContact::DIRICHLET))
          {
            double valu = (contact->get_boundary_value(POTENTIAL)
                + contact->get_inner_voltage()) / phi0;

            double valn = (contact->get_boundary_value(FERMIE)
                + contact->get_inner_voltage()) / phi0;

            double valp = (contact->get_boundary_value(FERMIH)
                + contact->get_inner_voltage()) / phi0;

            for (size_t n = 0; n < elem->n_nodes(); n++)
              if (elem->is_node_on_side(n, s))
              {
                if ((coupling & POISSON) &&
                    (contact->get_type(POTENTIAL) ==
                     ElectricalContact::DIRICHLET))
                  Ke.condense(n, n, -valu, Fe);

                if ((coupling & ECURRENT) &&
                    (contact->get_type(FERMIE) ==
                     ElectricalContact::DIRICHLET))
                  Ke.condense(n + n_dofs, n + n_dofs, -valn, Fe);

                if ((coupling & HCURRENT) &&
                    (contact->get_type(FERMIH) ==
                     ElectricalContact::DIRICHLET))
                  Ke.condense(n + 2 * n_dofs, n + 2 * n_dofs, -valp, Fe);
              }
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
    //jacobian->print_matlab("J.m");
  }
  else
  {
    residual->close();
    //residual->print_matlab("F.m");
    //ostringstream os;
    //os << "_" << __private_counter;
    //write_nodal_vector("residual" + os.str(), *residual);
    //write_nodal_vector("x" + os.str(), x);
    //__private_counter++;
  }


  perf_log.stop_event("assembly");
}



void
DriftDiffusion::save_data(const string& file)
{
  //using namespace boost::iostreams;
  //using namespace boost::iostreams::gzip;

  const MeshBase& mesh = get_mesh();
  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = static_cast<TiberNonlinearSystem&>(
      eq_sys.get_system(get_equation_system_name()));

  unsigned int sys_num = system.number();

  const NumericVector<Number>& solution = get_solution_vector();

  double phi0 = get_scaling().get_potential_scaling();

  ogzstream of(file.c_str());
  //ofstream of(file.c_str(), ios_base::out | ios_base::binary);
  //filtering_streambuf<output> out;
  //out.push(gzip_compressor());
  //out.push(of);

  // write contact voltages
  of << "<contacts>" << endl;

  ContactData sim_voltages(_boundary_currents);
  ContactData::iterator ctit(sim_voltages.begin());
  const ContactData::iterator ctend(sim_voltages.end());
  for ( ; ctit != ctend; ++ctit)
  {
    const Boundary* bd = ctit->first;
    // It's save to static_cast because we know there has to be an
    // ElectricalContact object
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(bd->get_boundary_properties(get_id()));
    of << bd->get_name() << " " << cnt->get_simulation_voltage() << endl;
  }

  of << "</contacts>" << endl;


  // write all variables
  of << "<variables>" << endl;
  Variable::iterator vit(Variable::begin());
  const Variable::iterator vend(Variable::end());
  for ( ; vit != vend; ++vit)
  {
    of << (*vit)->get_name() << " " << (*vit)->get_value_string() << endl;
  }
  of << "</variables>" << endl;

  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  unsigned int en_var = system.variable_number("fermi_e");
  unsigned int ep_var = system.variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  of << "<data>" << endl;

  Mesh::const_node_iterator it(mesh.active_nodes_begin());
  const Mesh::const_node_iterator end(mesh.active_nodes_end());
  for ( ; it != end; ++it)
  {
    const Node& node = *(*it);
    unsigned int dof_u  = node.dof_number(sys_num, u_var, 0);
    unsigned int dof_en = node.dof_number(sys_num, en_var, 0);
    unsigned int dof_ep = node.dof_number(sys_num, ep_var, 0);

    of << phi0 * solution(dof_u) << " " << phi0 * solution(dof_en) << " "
      << phi0 * solution(dof_ep) << endl << flush;
  }

  of << "</data>" << endl;
}



void
DriftDiffusion::load_data(const string& file)
{
  const MeshBase& mesh = get_mesh();
  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = static_cast<TiberNonlinearSystem&>(
      eq_sys.get_system(get_equation_system_name()));

  unsigned int sys_num = system.number();

  NumericVector<Number>& solution = get_solution_vector();

  double phi0 = get_scaling().get_potential_scaling();
  if (!get_option("scale_after_load", true))
    phi0 = 1.0;

  //ifstream is(file.c_str());
  igzstream is(file.c_str());
  if (!is.good()) throw InitFailedException("Bad datafile");

  string keyword("<contacts>");
  const streamsize bufsize = 256;
  char buf[bufsize];

  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(buf) != 0))
  {
    is.getline(buf, bufsize);
  }

  if (!is.good()) throw InitFailedException("Bad datafile");

  map<string, double> values;

  keyword = "</contacts>";
  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(buf) != 0))
  {
    istringstream ss(buf);
    string name;
    double value;
    ss >> name >> value;
    values[name] = value;
    is.getline(buf, bufsize);
  }


  ContactData sim_voltages(_boundary_currents);
  ContactData::iterator ctit(sim_voltages.begin());
  const ContactData::iterator ctend(sim_voltages.end());
  for ( ; ctit != ctend; ++ctit)
  {
    const Boundary* bd = ctit->first;
    // It's save to static_cast because we know there has to be an
    // ElectricalContact object
    ElectricalContact* cnt =
      static_cast<ElectricalContact*>(bd->get_boundary_properties(get_id()));

    cnt->set_simulation_voltage(values[bd->get_name()]);
    sim_voltages[bd] = cnt->get_simulation_voltage();
  }

  if (!is.good()) throw InitFailedException("Bad datafile");

  keyword = "<variables>";
  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(buf) != 0))
  {
    is.getline(buf, bufsize);
  }

  if (!is.good()) throw InitFailedException("Bad datafile");

  values.clear();

  keyword = "</variables>";
  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(buf) != 0))
  {
    istringstream ss(buf);
    string name;
    double value;
    ss >> name >> value;
    values[name] = value;
    is.getline(buf, bufsize);
  }

  map<string, double>::iterator vit(values.begin());
  const map<string, double>::iterator vend(values.end());
  for ( ; vit != vend; ++vit)
  {
    Variable::set_variable_value(vit->first, vit->second);
  }

  if (!is.good()) throw InitFailedException("Bad datafile");
  values.clear();

  keyword = "<data>";
  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(buf) != 0))
  {
    is.getline(buf, bufsize);
  }

  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_e");
  const unsigned int ep_var = system.variable_number("fermi_h");

  MeshBase::const_node_iterator it(mesh.active_nodes_begin());
  const MeshBase::const_node_iterator end(mesh.active_nodes_end());
  for ( ; it != end; ++it)
  {
    const Node& node = *(*it);
    unsigned int dof_u  = node.dof_number(sys_num, u_var, 0);
    unsigned int dof_en = node.dof_number(sys_num, en_var, 0);
    unsigned int dof_ep = node.dof_number(sys_num, ep_var, 0);

    if (!is.good()) throw InitFailedException("Bad datafile");

    is.getline(buf, bufsize);
    istringstream ss(buf);

    double u, en, ep;
    ss >> u >> en >> ep;

    solution.set(dof_u, u / phi0);
    solution.set(dof_en, en / phi0);
    solution.set(dof_ep, ep / phi0);
  }

  equilibrium_done(true);
  is_solved(true);
}



void
DriftDiffusion::write_nodal_vector(const string& filename, const NumericVector<double>& vec)
{

  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int nn  = mesh.n_nodes();

  vector<double> results(3 * nn);

  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

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

    ID subdomain = elem->subdomain_id();

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      unsigned int id = 3 * elem->node(n);
      results[id] = vec(dof_indices_u[n]);
      results[id + 1] = vec(dof_indices_en[n]);
      results[id + 2] = vec(dof_indices_ep[n]);
    }

  }

  DataOutput data_output(get_mesh(), "vtk");
  data_output.set_output_directory(get_output_directory());
  vector<string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  data_output.write_nodal_data(filename, results, names);

}

