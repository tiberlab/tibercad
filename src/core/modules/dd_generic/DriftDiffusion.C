// $Id: DriftDiffusion.C 4260 2016-05-27 11:23:40Z maufder $

// module includes
#include "DriftDiffusion.h"
#include "SimulationEnvironment.h"
#include "Scaling.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "Boundary.h"
#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"
#include "DDBulkModel.h"
#include "DDInterfaceModel.h"
#include "Constants.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSystem.h"
#include "SolveFailedException.h"
#include "Variable.h"
#include "FowlerNordheim.h"


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
#include "libmesh_logging.h"
#include "perf_log.h"
//#include "libMeshDefs.h"

#include "DataOutput.h"
#include "Messages.h"

// C++ includes
#include <fstream>

#include "TiberModule.h"

//
// Module interface
//


namespace
{
  int __private_counter;
}


using namespace std;
using namespace DriftDiffusionDefs;
using namespace libMesh;

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
    exact_newton(true),
    local_neutrality(false)
{
}






DriftDiffusion::DriftDiffusion(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true),
    _useparticle('b'),
    _reference_potential(0.0),
    _rstf(NULL)
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


  DDBulkModel* model =
    DDBulkModel::create(modelname, mat, options);

  if (model == NULL)
    throw ModelErrorException(
        "DriftDiffusion: No such physical model: " + modelname);

  return model;
}



PhysicalModel*
DriftDiffusion::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{

  PhysicalModel* model = NULL;

  model = DDInterfaceModel::create(boundary, options);

  return model;
}






void
DriftDiffusion::compute_scaling(Scaling::ScalingType type)
{

  // we calculate in cm!
  double mesh_units = 100 * get_mesh_units();
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
  MeshBase::const_element_iterator el = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    assert(_device->get_material(elem->subdomain_id()) != NULL);
    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

    sc->set_coordinates(elem->centroid());
    sc->set_potentials(sc->get_equilibrium_fermi_level());
    sc->set_electric_field(libMesh::RealGradient(0));
    for (auto&& cp : sc->get_carrier_properties())
      sc->set_grad_fermi(cp.first, libMesh::RealGradient(0));
    sc->reinit(elem);

    sc->calculate_densities();
    sc->calculate_traps();
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

    const libMesh::RealTensor& eps_tens = sc->get_relative_permittivity();
    eps0 = (eps0 > eps_tens(0,0)) ? eps0 : eps_tens(0,0);

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

      pair<Point, Point> bbox(get_environment().get_bounding_box());
      Point dia(bbox.second - bbox.first);

      x0 = (dia(0) > dia(1)) ? dia(0) : dia(1);
      x0 = (x0 > dia(2)) ? x0 : dia(2);

      break;
  }

  get_scaling().set_scaling_type(type);
  get_scaling().set_potential_scaling(phi0);
  get_scaling().set_length_scaling(x0 * mesh_units);
  vector<double> values = {mu0, C0};
  this->get_communicator().max(values);
  mu0 = values[0];
  C0  = values[1];
  get_scaling().set_mobility_scaling(mu0 > 0 ? mu0 : 1.0);
  get_scaling().set_density_scaling(C0 > 0 ? C0 : 1.);
}




void
DriftDiffusion::set_electron_fermi_level(double Ef_n)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = get_solution_vector();

  const unsigned int var = system.variable_number("fermi_e");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = Ef_n / phi0;

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = this->active_local_elements_begin();
  const MeshBase::element_iterator end = this->active_local_elements_end();

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
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = get_solution_vector();

  const unsigned int var = system.variable_number("fermi_h");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = Ef_p / phi0;

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = this->active_local_elements_begin();
  const MeshBase::element_iterator end = this->active_local_elements_end();

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
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = get_solution_vector();

  const unsigned int var = system.variable_number("potential");
  const double phi0 = get_scaling().get_potential_scaling();
  double level = -pot / phi0;

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = this->active_local_elements_begin();
  const MeshBase::element_iterator end = this->active_local_elements_end();

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
DriftDiffusion::find_dielectric_boundary_nodes(void)
{
  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = this->active_local_elements_begin();
  const MeshBase::element_iterator end = this->active_local_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* el = *it;

    DDBulkModel* sc = get_bulk_model<DDBulkModel>(el);

    // we are only interested in boundaries between semiconductor/dielectric
    if (sc->is_dielectric())
    {
      for (unsigned s = 0; s < el->n_sides(); s++)
      {
        if (get_environment().is_inner_boundary(ElementSide(el, s)))
        {
          // get the model of the neighbor element
          DDBulkModel* scn = get_bulk_model<DDBulkModel>(el->neighbor(s));


          // if neighbor is not dielectric we record it
          if (!scn->is_dielectric())
          {
            libMesh::UniquePtr<Elem> side(el->build_side(s));
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
      ElectricalContact* bd = dynamic_cast<ElectricalContact*>(
          it->first->models_begin()->second);

      if (bd != NULL)
      {
        double voltage = bd->get_simulation_voltage();

        _voltages[it->first] = voltage;

        if (voltage != 0.0)
          equilibrium = false;
      }
    }
  }



  if (!equilibrium_done())
  {
    solve_equilibrium();

    build_local_scaling();

    // if we would repeat the equilibrium simulation, we can stop now
    if (equilibrium)
      return;
  }

  // set the old solution
  //EquationSystems& es = get_equation_systems();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();
  get_solution_vector().close();
  system.get_vector("old_sol") = get_solution_vector();


  int coupling = get_my_options().coupling;

/***
  if (get_options().find_option("guess_el_qfermi") ||
      get_options().find_option("guess_hl_qfermi"))
  {
    get_my_options().coupling = POISSON;

    if (get_options().find_option("guess_el_qfermi"))
      set_electron_fermi_level(_el_qfermi_guess);

    if (get_options().find_option("guess_hl_qfermi"))
      set_hole_fermi_level(_hl_qfermi_guess);

    do_newton();

    get_my_options().coupling = coupling;
  }
***/

  //set_dirichlet_bc();
  if (get_option("use_weight", false))
    calculate_weights();

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
  catch (::SolverException& e)
  {
    string msg = "solve failed (" +
        string(e.what()) + ")";
    throw SolveFailedException(msg);
  }

  // NOTE we calculate the local scaling factors AFTER, because otherwise
  // new results of other coupled models (thermal, quantum) may disturb
  // the calculation
  build_local_scaling();

  get_my_options().coupling = coupling;

  // calculate the currents to print them on screen
  calculate_currents();
  calculate_iqe();
  calculate_surface_recombination();


  ContactData::iterator it(_boundary_currents.begin());
  const ContactData::iterator end(_boundary_currents.end());

  int width = 20;
  {
    Messages::newline();
    ostringstream os;
    os << "contact name:";
    os.width(width - os.tellp());
    os << "";
    os << "voltage:";
    os.width(2 * width - os.tellp());
    os << "";
    os << "current:";
    os << endl;
    Messages::info(os.str());
  }

  for (it = _boundary_currents.begin(); it != end; ++it)
  {
    DDInterfaceModel* ifmod =
        static_cast<DDInterfaceModel*>(it->first->models_begin()->second);

    // we print only contacts with a current
    if (!ifmod->has_current())
      continue;

    ostringstream os;
    os << setprecision(6);
    os << it->first->get_name();
    ElectricalContact* cnt = dynamic_cast<ElectricalContact*>(ifmod);
    if (cnt != NULL)
    {
      if (cnt->get_simulation_voltage() < 0)
        os.width(width - os.tellp() - 1);
      else
        os.width(width - os.tellp());
      os << "";
      os << cnt->get_simulation_voltage();
    }
    else
    {
      os.width(width - os.tellp());
      os << "" << "-";
    }


    if (it->second < 0)
      os.width(2 * width - os.tellp() - 1);
    else
      os.width(2 * width - os.tellp());
    os << "";
    os << it->second * it->first->get_area_factor();

    Messages::info(os.str());
  }

  if (plot_solution(IQE))
  {
    ostringstream os;
    os << setprecision(6);
    os << "IQE: " << _iqe;

    Messages::newline();
    Messages::info(os.str());

  }

  calculate_field_emission();

  // calculate the mean fermi levels on the contacts
  // TODO this one has some problem in parallel
  //calculate_mean_fermi_levels();

}


void
DriftDiffusion::calculate_weights(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = get_solution_vector();
  libMesh::NumericVector<Number>& oldsol = system.get_vector("old_sol");
  libMesh::NumericVector<Number>& weight = system.get_vector("weight");

  const unsigned int var_u = system.variable_number("potential");
  // fermi potential variable numbers are defined within element loop
  // since (in future) different variables can ben set in different regions

  const double phi0 = get_scaling().get_potential_scaling();
  const double C0 = get_scaling().get_density_scaling();

  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = this->active_local_elements_begin();
  const MeshBase::element_iterator end = this->active_local_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

    assert(sc != NULL);
    sc->reinit(elem);

    //Get variables for fermi potentials
    map<unsigned int, string> var_q;
    for (auto&& cp : sc->get_carrier_properties())
      var_q.insert( make_pair(system.variable_number(cp.first), cp.first) );

    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      sc->set_coordinates(elem->point(i));

      unsigned int dofu =
        elem->get_node(i)->dof_number(system.number(), var_u, 0);
      map<unsigned int, unsigned int> dofq;
      for (auto var : var_q)
        dofq.insert( make_pair(var.first, 
                               elem->get_node(i)->
                                       dof_number(system.number(), var.first, 0)) );

      double u = solution(dofu);
      double oldu = oldsol(dofu);
      map<unsigned int, double> q, oldq;
      for (auto& dof : dofq)
      {
        q.insert( make_pair(dof.first, solution(dof.second)) );
        oldq.insert( make_pair(dof.first, oldsol(dof.second)) );
      }

      sc->set_el_potential(phi0 * u);
      sc->set_old_el_potential(phi0 * oldu);
      for (auto var : var_q)
      {
        sc->set_fermi_potential(var.second, phi0 * q[var.first]);
        sc->set_old_fermi_potential(var.second, phi0 * q[var.first]);
      }

      sc->calculate_densities();

      weight.set(dofu, 1);

      for (auto var : var_q)
      {
        double dens = min(1.0, log10(1 + sc->get_q_density(var.second)) / 18);// < 1e7 ? 0 : 1;
        weight.set(dofq[var.first], dens);
      }

    }
  }
  weight.close();
  system.set_weight(&weight, TiberEqSystem::l2_NORM);
  system.set_weight(&weight, TiberEqSystem::MAX_NORM);
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

  // update the number of components of the band edge variables
  {
    size_t num_eb = 0;
    size_t num_hb = 0;

    const set<PhysicalModel*>& pm = get_physical_models();
    set<PhysicalModel*>::const_iterator it(pm.begin());
    set<PhysicalModel*>::const_iterator end(pm.end());
    for ( ; it != end; ++it)
    {
      DDBulkModel* sc =
          static_cast<DDBulkModel*>(*it);

      // carrier bands
      map<string, vector<double>> eb;
      map<string, vector<double>> hb;
      for (auto&& cp: sc->get_carrier_properties())
      {
        vector<double> bands;
        cp.second->get_bands(bands);
        if (cp.second->get_carrier_type() == 'e')
        {
          eb.insert(make_pair(cp.first, bands));
          num_eb = max(num_eb, eb.size());
        }
        else
        {
          hb.insert(make_pair(cp.first, bands));
          num_hb = max(num_hb, hb.size());
        }
      }

    }

    declare_solution(ElectronBands, NTUPLE, CELL, "eV", num_eb);
    declare_solution(HoleBands, NTUPLE, CELL, "eV", num_hb);
  }


  // first we have to compute the scaling
  compute_scaling(get_my_options().scaling_type);


  ModelOptions& solveropts = get_solver_options();
  int max_it = solveropts.get_option("max_iterations", 15);
  if (max_it < 100)
  {
    Messages::info("Setting max_iterations for nonlinear solver to 150");
    solveropts.set_option("max_iterations", 150);
  }

  int coupling = get_my_options().coupling;
  get_my_options().coupling = POISSON;

  // backup the simulation voltages and set all to zero
  ContactData sim_voltages(_boundary_currents);
  ContactData::iterator it(sim_voltages.begin());
  const ContactData::iterator end(sim_voltages.end());
  for ( ; it != end; ++it)
  {
    const Boundary* bd = it->first;
    vector<ID> ids;
    bd->get_region_ids(ids);
    for (unsigned int i = 0; i < ids.size(); i++)
    {
      // it must exist
      MaterialBoundary* mb = get_environment().get_device().get_boundary_object(ids[i]);
      ElectricalContact* cnt =
          dynamic_cast<ElectricalContact*>(mb->get_model(get_id()));

      if (cnt != NULL)
      {
        sim_voltages[bd] = cnt->get_simulation_voltage();
        cnt->set_simulation_voltage(0.0);
      }
    }
  }

  // make a rough guess
  guess_equilibrium();

  build_local_scaling();



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
    vector<ID> ids;
    bd->get_region_ids(ids);
    for (unsigned int i = 0; i < ids.size(); i++)
    {
      // it must exist
      MaterialBoundary* mb = get_environment().get_device().get_boundary_object(ids[i]);
      ElectricalContact* cnt =
          dynamic_cast<ElectricalContact*>(mb->get_model(get_id()));

      if (cnt != NULL)
        cnt->set_simulation_voltage(sim_voltages[bd]);
    }
  }

  // reset the coupling
  get_my_options().coupling = coupling;

  solveropts.set_option("max_iterations", max_it);


  // compute the reference potential
  compute_reference_potential();

}



void
DriftDiffusion::compute_reference_potential(void)
{
  if (get_my_options().reference_contact != "")
  {
    SimulationEnvironment& si = get_environment();
    std::set<const Node*> nodelist;
    si.get_boundary_nodes(get_my_options().reference_contact, nodelist);
    if (nodelist.size() > 0)
    {
      TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

      const libMesh::NumericVector<Number>& solution = get_solution_vector();
      const unsigned int system_number = system.number();
      const unsigned int u_var = system.variable_number("potential");

      const Node* node = *nodelist.begin();
      const unsigned int n_dof = node->dof_number(system_number, u_var, 0);
      _reference_potential = solution(n_dof) * get_scaling().get_potential_scaling();
    }
  }

}



void
DriftDiffusion::calculate_iqe(void)
{
  _iqe = 0;

  set<ID> active_regs = get_region_ids();
  bool recomb_only = false;

  ModelOptions::submodel_iterator opts(get_options().submodels_begin("iqe"));
  if (opts != get_options().submodels_end("iqe"))
  {
    recomb_only = (opts->second).get_option("use_recombinations_only", recomb_only);
    string reg = (opts->second).get_option("active_regions", "all");
    set<ID> tmp;
    get_environment().get_device().extract_physical_regions(reg, tmp);

    set<ID>::iterator it(active_regs.begin());
    set<ID>::iterator end(active_regs.end());
    while (it != end)
    {
      set<ID>::iterator tmpit(it);
      ++it;
      if (!tmp.count(*tmpit)) active_regs.erase(tmpit);
    }
  }

  ID rec_id = get_solution_id("eDirectRecombination");
  ID netrec_id = get_solution_id("eNetRecombination");
  ID srh_id = get_solution_id("eSRHRecombination");
  ID aug_id = get_solution_id("eAugerRecombination");
  double Rtot = 0.0, Rsrh = 0.0, Raug = 0.0;

  if (rec_id == INVALID_ID)
  {
    if (plot_solution(IQE))
    {
      Messages::warning("Cannot calculate IQE: no direct recombination model present.");
    }
    return;
    //throw RuntimeException("Recombination model \'DirectRecombination\' "
    //    "is not defined in Drift-Diffusion simulation.");
  }

  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = system.variable_number("potential");

  libMesh::FEType fe_type = system.variable_type(u_var);

  const Options& params = get_my_options();

  libMeshEnums::Order integration_order = params.integration_order;

  // the finite element
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
  libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(
        params.quadrature_type, dim, integration_order));
  fe->attach_quadrature_rule(qrule.get());

  const vector<Real>& JxW = fe->get_JxW();

  // the finite element for boundary integration
  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));

  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;

  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
        params.quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());

  // references to boundary-specific data that will be used to
  // assemble the system.
  // Data will be given for each quadrature point.
  //
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  //
  const vector<vector<libMesh::RealGradient> >&  dphi_face = fe_face->get_dphi();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point_face = fe_face->get_xyz();
  //
  const vector<Point>& face_normals = fe_face->get_normals();
  //
  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW_face = fe_face->get_JxW();


  MeshBase::const_element_iterator el =
      this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();

    if (!active_regs.count(subdomain)) continue;

    fe->reinit(elem);

    map<ID, vector<double> > datamap;
    datamap[rec_id] = vector<double>(qrule->n_points());
    datamap[netrec_id] = vector<double>(qrule->n_points());
    datamap[srh_id] = vector<double>(qrule->n_points());
    datamap[aug_id] = vector<double>(qrule->n_points());
    vector<double>& data = datamap[rec_id];
    vector<double>& nrdata = datamap[netrec_id];
    vector<double>& srhdata = datamap[srh_id];
    vector<double>& augdata = datamap[aug_id];

    DriftDiffusion::get_solution_secure(elem, datamap, qrule->get_points());

    double iqe_el = 0;

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {
      iqe_el += JxW[qp] * data[qp];
      Rtot += JxW[qp] * nrdata[qp];
      Rsrh += JxW[qp] * srhdata[qp];
      Raug += JxW[qp] * augdata[qp];
    }

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {

      Material* mat = get_material(elem);
      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      //bool true_boundary = environment.is_outer_boundary(ElementSide(elem, s));

      //if ((sm != NULL) || true_boundary)
      //{
      //  fe_face->reinit(elem, s);

      //  int phi_size = phi_face.size();

        // now integrate to include von Neumann and mixed type BCs
        // and polarization
      //  for (unsigned int qp = 0; qp < qface->n_points(); qp++)
      //  {
      //  }
      //}
    }

    _iqe += iqe_el;
  }

  vector<double> recs = {_iqe, Rsrh, Raug, Rtot};
  this->get_communicator().sum(recs);
  _iqe = recs[0];
  Rsrh = recs[1];
  Raug = recs[2];
  Rtot = recs[3];

  // 2015-03-02: we take now only direct, SRH and Auger for the definition
  // this eliminates problems in presence of optical generation
  Rtot = _iqe + Rsrh + Raug;

  ostringstream rec;
  rec << "Rrad = " << _iqe * Constants::e << 
    " Rsrh = " << Rsrh * Constants::e <<
    " Raug = " << Raug * Constants::e << 
    " Rtot = " << Rtot * Constants::e << "\n";
  Messages::info(rec.str());

  if (recomb_only)
  {
    _iqe /= Rtot;
  }
  else
  {
    double current = 0;

    // now take the total outflowing current
    ContactData::const_iterator it(_boundary_currents.begin());
    for ( ; it != _boundary_currents.end(); ++it)
    {
      DDInterfaceModel* ifmod =
          static_cast<DDInterfaceModel*>(it->first->models_begin()->second);

      // we print only contacts with a current
      if (!ifmod->has_current())
        continue;

      if (it->second > 0)
        current += it->second;
    }

    if (current > 0)
      _iqe *= Constants::e / current;
    else
      _iqe = 0;

  }
}


void
DriftDiffusion::guess_equilibrium(void)
{

  // equation system needs to be active
  rebuild_equation_system();

  TiberNonlinearSystem& poisson = get_equation_system<TiberNonlinearSystem>();

  const unsigned int u_var = poisson.variable_number("potential");
  const libMesh::DofMap& dof_map_u = poisson.get_dof_map();
  vector<unsigned int> dof_indices_u;

  libMesh::NumericVector<Number>& solution_u = poisson.get_solution_vector();
  solution_u.close();
  //solution_u.zero();

  MeshBase::const_element_iterator el =
                                  this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  this->active_local_elements_end();

  const double phi0 = get_scaling().get_potential_scaling();

  const unsigned int nn  = get_mesh().n_nodes();
  vector<unsigned short int> node_conn(nn);
  // Get the number of elements that share each node.  We will
  // compute the average value at each node.
  {
    MeshBase::const_element_iterator it =
      this->active_local_elements_begin();
    const MeshBase::const_element_iterator end =
      this->active_local_elements_end();

    for ( ; it != end; ++it)
      for (unsigned int n = 0; n < (*it)->n_nodes(); n++)
	node_conn[(*it)->node(n)]++;

    this->get_communicator().sum(node_conn);
  }

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();
    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

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

  get_environment().prepare_for_solve();
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

  if (_do_local_scaling)
    os << "using local scaling";

  Messages::info(os.str());

}





void
DriftDiffusion::parse_const_options(void)
{

  const ModelOptions& opts = get_options();
  Options& myopts = get_my_options();

  string method =  opts.get_option("current_integration_method", "rstf");
  if (method == "rstf")
    myopts.current_calculation = RSTF;
  else if (method == "compact_rstf")
    myopts.current_calculation = RSTF_COMPACT;
  else if (method == "surface_integral")
    myopts.current_calculation = SURFINT;
  else throw InitFailedException("Unknown current integration method: "
      + method + " in DriftDiffusion.");

  string scaling = opts.get_option("scaling", "");
  if (scaling == "demari")
    myopts.scaling_type = Scaling::DEMARI;
  else if (scaling == "none")
    myopts.scaling_type = Scaling::NONE;
  else
    myopts.scaling_type = Scaling::UNITS;

  _do_local_scaling = opts.get_option("local_scaling", true);

  //string qrule = get_mesh().mesh_dimension() == 1 ? "trapez" : "gauss";
  string qrule = "trapez";
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
    if (opts.get_option("local_equilibrium", true))
      _useparticle = 'e';
  }
  else if (coupling == "holes")
  {
    myopts.coupling = HCURRENT | POISSON;
    if (opts.get_option("local_equilibrium", true))
      _useparticle = 'h';
  }
  else if (coupling == "current")
    myopts.coupling = CURRENTS;

  // for Poisson only simulations one may want to enforce local
  // charge neutrality
  myopts.local_neutrality = get_option("enforce_local_charge_neutrality",
      myopts.local_neutrality);
  if (myopts.local_neutrality && (myopts.coupling != POISSON))
  {
    Messages::warning("You try to enforce local charge neutrality in " +
        get_name() + " even if you solve for continuity equations.");
  }


  myopts.mesh_refinement = opts.get_option("mesh_refinement",
      myopts.mesh_refinement);

  myopts.exact_newton = opts.get_option("exact_newton", myopts.exact_newton);

  get_parameter("guess_el_qfermi", _el_qfermi_guess);
  get_parameter("guess_hl_qfermi", _hl_qfermi_guess);

  myopts.reference_contact = opts.get_option("reference_contact", "");

  string s(opts.get_option("default_boundary_condition", "zero_field"));
  if (s == "zero_field")
  {
    myopts.default_boundary_condition = ZEROFIELD;
  }
  else if (s == "zero_displacement")
  {
    myopts.default_boundary_condition = ZERODISPLACEMENT;
  }
  else
    Messages::warning(s + " is unknown Poisson boundary condition.");

}





void
DriftDiffusion::rebuild_equation_system(void)
{
  if (!_rebuild_eq_system) return;


  //EquationSystems& equation_systems = get_equation_systems();
  clear_systems();

  ModelOptions::submodel_iterator linit(
      get_solver_options().submodels_begin("linear_solver"));

  if (linit == get_solver_options().submodels_end("linear_solver"))
  {
    get_solver_options().add_submodel("linear_solver", ModelOptions());
    linit = get_solver_options().submodels_begin("linear_solver");
  }

  ModelOptions& linopts = linit->second;

  if ((linopts.get_name() == "") ||
      (linopts.get_name() == "petsc"))
  {
    // default is bcgsl
    if (!linopts.find_option("method"))
      linopts["method"] = "bcgsl";

    // in 1D bcgs seems to work better than bcgsl
    const unsigned int dim = get_mesh().mesh_dimension();
    if ((dim == 1) && (linopts["method"] == "bcgsl"))
      linopts["method"] = "bcgs";

    if (!linopts.find_option("preconditioner"))
    {
      if (dim < 3)
        linopts["preconditioner"] = "lu";
      else
        linopts["preconditioner"] = "ilu";
    }

    if (linopts.get_option("absolute_tolerance", -1.0) < 0)
      linopts["absolute_tolerance"] = "1e-15";
  }



  ModelOptions& solveropts = get_solver_options();
  if (solveropts.get_option("absolute_tolerance", -1.0) < 0)
    solveropts["absolute_tolerance"] = "1e-15";


  // the coupled DD system
  create_equation_system("nonlinear");
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  system.attach_assembly_routine(assemble_system);

  ModelOptions::submodel_iterator physit (get_options().submodels_begin("Physics"));
  ModelOptions& physopts = physit->second;

  ModelOptions::submodel_iterator itc (physopts.submodels_begin("carrier"));
  ModelOptions::submodel_iterator end_itc (physopts.submodels_end("carrier"));

  system.add_variable("potential", libMeshEnums::FIRST);
  _variables.resize(0);
  for ( ; itc != end_itc; ++itc)
  {
    string name = (itc->second).get_option("name", "");
    _variables.push_back(name);
    system.add_variable(name, libMeshEnums::FIRST);
  }

  system.add_vector("old_sol");
  system.add_vector("weight");
  system.add_vector("scaling");
  //system.add_matrix("Preconditioner");



  // finally initialize the newly created system
  system.init();


  _rebuild_eq_system = false;

}


DriftDiffusion::RSTFSys*
DriftDiffusion::RSTFSys::create(DriftDiffusion* dd)
{

  libMesh::EquationSystems& es = dd->get_equation_systems();
  DriftDiffusion::RSTFSys* sys = NULL;
  sys = &(es.add_system<RSTFSys>("__DD_rstf"));
  if (sys == NULL)
    throw InitFailedException("Fatal error in DriftDiffusion. "
        "Cannot create RSTF system");

  sys->_dd = dd;

  sys->add_variable("u", libMeshEnums::FIRST, &dd->get_region_ids());
  // for each contact, we add a vector for the Ramo-Shockley test function
  ContactData::iterator it = dd->_boundary_currents.begin();
  for (int i = 0 ; it != dd->_boundary_currents.end(); ++it, ++i)
  {
    ostringstream os;
    os << "rstf" << i;
    sys->add_vector(os.str(), true, GHOSTED);
    sys->_boundaries[it->first] = i;
  }

  sys->assemble_before_solve = false;

  return sys;
}


void
DriftDiffusion::RSTFSys::solve(void)
{
  assemble();

  // we need to solve for every contact
  map<const Boundary*, int>::iterator bdit(_boundaries.begin());
  for ( ; bdit != _boundaries.end(); ++bdit)
  {
    int i = bdit->second;

    ostringstream os;
    os << "rstf" << i;

    get_vector(os.str()).close();
    *rhs = get_vector(os.str());

    LinearImplicitSystem::solve();

    get_vector(os.str()) = *current_local_solution;

  }

  //plot();
}


libMesh::NumericVector<double>*
DriftDiffusion::RSTFSys::get_testfunction(int i)
{
  ostringstream os;
  os << "rstf" << i;
  return &get_vector(os.str());
}

libMesh::NumericVector<double>*
DriftDiffusion::RSTFSys::get_testfunction(const Boundary* bd)
{
  libMesh::NumericVector<double>* vec = NULL;

  map<const Boundary*, int>::iterator bdit(_boundaries.find(bd));
  if (bdit != _boundaries.end())
  {
    ostringstream os;
    os << "rstf" << bdit->second;
    vec = &get_vector(os.str());
    vec->close();
  }

  return vec;
}


void
DriftDiffusion::RSTFSys::user_assembly(void)
{
  const unsigned int var = variable_number("u");

  vector<libMesh::NumericVector<double>*> rhsides(_boundaries.size());
  map<const Boundary*, int>::iterator bdit(_boundaries.begin());
  for ( ; bdit != _boundaries.end(); ++bdit)
  {
    int i = bdit->second;
    ostringstream os;
    os << "rstf" << i;
    rhsides[i] = &get_vector(os.str());
  }

  libMesh::DofMap& dof_map = get_dof_map();

  libMesh::FEType fe_type = variable_type(var);

  const MeshBase& mesh = get_mesh();
  unsigned int dim = mesh.mesh_dimension();

  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
  libMesh::QGauss qrule(dim, SECOND);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  std::vector<unsigned int> dof_indices;

  libMesh::DenseMatrix<Number> Ke;
  vector<libMesh::DenseVector<Number>> Fe(_boundaries.size());

  const double penalty = 1e6;

  MeshBase::const_element_iterator el(_dd->active_local_elements_begin());
  const MeshBase::const_element_iterator end_el(_dd->active_local_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs   = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (int i = 0; i < _boundaries.size(); ++i)
      Fe[i].resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
      for (unsigned int i = 0; i < n_dofs; i++)
        for (unsigned int j = 0; j < n_dofs; j++)
          Ke(i, j) += JxW[qp] * dphi[i][qp] * dphi[j][qp];


    // loop over the sides for boundary conditions
    // NOTE we use penalty-method here
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      Boundary* bd = _dd->get_environment().get_boundary(ElementSide(elem, s));

      if (bd != NULL)
      {
        map<const Boundary*, int>::iterator bdit(_boundaries.find(bd));
        if (bdit != _boundaries.end())
        {
          // the internal number of the boundary
          int bdid = bdit->second;

          for (unsigned int i = 0; i < elem->n_nodes(); i++)
          {
            //if (elem->is_node_on_side(i, s))
            {
              Ke(i, i) += penalty;
              Fe[bdid](i) += penalty;
            }
          }
        }
      }
    }

    // apply dof constraints
    // TODO needs to be split and uncommented for mesh refinement
    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    matrix->add_matrix(Ke, dof_indices);
    map<const Boundary*, int>::iterator bdit(_boundaries.begin());
    for ( ; bdit != _boundaries.end(); ++bdit)
    {
      rhsides[bdit->second]->add_vector(Fe[bdit->second], dof_indices);
    }

  }
}



void
DriftDiffusion::RSTFSys::plot(void)
{
  // we return immediately if nothing is to be printed
  //if (!_do_plot) return;

  const MeshBase& mesh = get_mesh();

  DataOutput data_output(mesh, TiberCad::get_output_format());
  data_output.set_output_directory(TiberCad::get_output_dir());


  vector<double> sol;
  vector<string> solname;

  build_nodal_results(sol, solname);
  data_output.write_nodal_data("__DD_rstf", sol, solname);

}


void
DriftDiffusion::RSTFSys::build_nodal_results(vector<double>& results,
    vector<string>& legend)
{
  legend.resize(_boundaries.size());

  vector<libMesh::NumericVector<double>*> rhsides(_boundaries.size());
  map<const Boundary*, int>::iterator bdit(_boundaries.begin());
  for ( ; bdit != _boundaries.end(); ++bdit)
  {
    int i = bdit->second;
    ostringstream os;
    os << "rstf" << i;
    legend[i] = os.str();
    rhsides[i] = &get_vector(os.str());
  }

  //const unsigned int s = number();
  //const unsigned int var = variable_number("u");
  const MeshBase& mesh = get_mesh();

  // TODO this iterators would not work on subdomains
  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  unsigned int number_of_points = 0;
  for ( ; nd != nd_end; ++nd)  number_of_points++;

  results.resize(number_of_points * _boundaries.size(), 0.0);

  MeshBase::const_element_iterator it(_dd->active_local_elements_begin());
  const MeshBase::const_element_iterator
    end(_dd->active_local_elements_end());

  vector<unsigned int> dof_indices;
  libMesh::DofMap& dof_map = get_dof_map();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      unsigned int id = elem->node(n) * _boundaries.size();
      map<const Boundary*, int>::iterator bdit(_boundaries.begin());
      for ( ; bdit != _boundaries.end(); ++bdit)
      {
        int i = bdit->second;
        results[id + i]  =  (*rhsides[i])(dof_indices[n]);
      }
    }
  }
}

void
DriftDiffusion::prepare_rstf(void)
{
  _rstf = RSTFSys::create(this);
  _rstf->init();
  _rstf->solve();

  _rstf->matrix->clear();
}




void
DriftDiffusion::do_init(void)
{

  _device = &get_environment().get_device();

  parse_const_options();

  find_dielectric_boundary_nodes();


  get_environment().update_boundary_element_map();


  // check the boundaries if they are internal or not
  const set<PhysicalModel*>& pm = get_interface_models();
  set<PhysicalModel*>::const_iterator it(pm.begin());
  const set<PhysicalModel*>::const_iterator end(pm.end());
  for ( ; it != end; ++it)
  {
    DDInterfaceModel* mod = static_cast<DDInterfaceModel*>(*it);

    const MaterialBoundary* bd =
        static_cast<const MaterialBoundary*>(mod->get_owner());


    // register the contact if it is a real contact (with current)
    if (mod->has_current())
    {
      const Boundary* bnd = get_environment().get_boundary(mod->get_name());
      if (bnd != NULL)
      {
        _boundary_currents[bnd] = 0.0;
        _voltages[bnd] = 0.0;
      }
    }


    // check if the boundary is internal or not
    ID idA = bd->get_id_A();
    ID idB = bd->get_id_B();

    if (get_environment().contains_region(idA) &&
        get_environment().contains_region(idB))
      mod->internal_bondary(true);
  }

  rebuild_equation_system();

}




void
DriftDiffusion::do_reinit(void)
{
  const set<PhysicalModel*>& pm = get_physical_models();
  set<PhysicalModel*>::const_iterator it(pm.begin());
  set<PhysicalModel*>::const_iterator end(pm.end());

  for ( ; it != end; ++it)
  {
    DDBulkModel* sc =
        static_cast<DDBulkModel*>(*it);

    switch (_useparticle)
    {
      case 'e':
        sc->set_coupling_type(DriftDiffusionDefs::ELECTRONS);
        break;

      case 'h':
        sc->set_coupling_type(DriftDiffusionDefs::HOLES);
        break;

      default:
        sc->set_coupling_type(DriftDiffusionDefs::BOTH);
        break;
    }
  }
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

  // the correct number of components will be inserted afterwards
  declare_solution(ElectronBands, NTUPLE, CELL, "eV", 1);
  declare_solution(HoleBands, NTUPLE, CELL, "eV", 1);
  if (plot_solution("BandEdges"))
  {
    add_plot_variable(ElectronBands);
    add_plot_variable(HoleBands);
  }

  declare_solution(Polarization, VECTOR, CELL, "C/m^2");

  declare_solution(eDensity, REAL, NODES, "cm^-3");
  declare_solution(hDensity, REAL, NODES, "cm^-3");
  if (plot_solution("Density"))
  {
    add_plot_variable(eDensity);
    add_plot_variable(hDensity);
  }

  declare_solution(eMobility, REAL, NODES, "cm^2/(V*s)");
  declare_solution(hMobility, REAL, NODES, "cm^2/(V*s)");
  if (plot_solution("Mobility"))
  {
    add_plot_variable(eMobility);
    add_plot_variable(hMobility);
  }

  declare_solution(eConductivity, REAL, NODES, "S/cm");
  declare_solution(hConductivity, REAL, NODES, "S/cm");
  if (plot_solution("Conductivity"))
  {
    add_plot_variable(eConductivity);
    add_plot_variable(hConductivity);
  }

  declare_solution(CurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(eCurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(hCurrentDensity, VECTOR, CELL, "A/cm^2");
  //declare_solution(eFlux, VECTOR, CELL, "1/(s*cm^2)");
  //declare_solution(hFlux, VECTOR, CELL, "1/(s*cm^2)");

  declare_solution(IonizedDonors, REAL, NODES, "cm^-3");
  declare_solution(IonizedAcceptors, REAL, NODES, "cm^-3");

  declare_solution(IonizedElectronTraps, REAL, NODES, "cm^-3");
  declare_solution(IonizedHoleTraps, REAL, NODES, "cm^-3");

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

  if (plot_solution("NetRecombination"))
  {
    add_plot_variable(eNetRecombination);
    add_plot_variable(hNetRecombination);
  }
  declare_solution(eNetRecombination, REAL, NODES, "1/(s*cm^3)");
  declare_solution(hNetRecombination, REAL, NODES, "1/(s*cm^3)");
  // add the single recombination rates
  {
    //size_t num_cb = 1;
    //size_t num_vb = 1;

    bool plot_erec = plot_solution(eNetRecombination);
    bool plot_hrec = plot_solution(hNetRecombination);

    const set<PhysicalModel*>& pm = get_physical_models();
    set<PhysicalModel*>::const_iterator it(pm.begin());
    set<PhysicalModel*>::const_iterator end(pm.end());

    for ( ; it != end; ++it)
    {
      DDBulkModel* sc =
          static_cast<DDBulkModel*>(*it);

      // the conduction bands
      //const vector<double>& cb = sc->get_conduction_bands();
      //num_cb = max(num_cb, cb.size());

      // the valence bands
      //const vector<double>& vb = sc->get_valence_bands();
      //num_vb = max(num_vb, vb.size());

      // we let it create the recombination models first
      sc->create_recombination_models();

      vector<ID> ids;
      int n = sc->get_net_recombination_rate_IDs(ids);

      for (int i = 0; i < n; i++)
      {
        ID eid = static_cast<ID>(eNetRecombination) + ids[i];
        ID hid = static_cast<ID>(hNetRecombination) + ids[i];
        const std::string& name =
            sc->get_recombination_model(ids[i])->get_default_name();
        // if recombination should be plotted, add it also to the plot variables
        if (plot_erec) add_plot_variable("e" + name);
        if (plot_hrec) add_plot_variable("h" + name);
        declare_solution_ext("e"+name, eid, SolutionDescriptor::REAL,
            SolutionDescriptor::NODES, "1/(s*cm^3)");
        declare_solution_ext("h"+name, hid, SolutionDescriptor::REAL,
            SolutionDescriptor::NODES, "1/(s*cm^3)");

        _recombination_ids.insert(ids[i]);
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

    //declare_solution(ContactCurrent, REAL, GLOBAL, units);
    //add_alias("ContactCurrents", ContactCurrent);
    bool plot_curr = plot_solution("ContactCurrent");
    plot_curr |= plot_solution("ContactCurrents");

    // we put them first into a set so they will be ordered
    // alphabetically
    set<string> cnt_set;

    SimulationEnvironment::BoundaryIterator it(get_environment().boundaries_begin());
    const SimulationEnvironment::BoundaryIterator end(get_environment().boundaries_end());

    for ( ; it != end; ++it)
    {
      const Boundary* bd = (*it);
      if (bd != NULL)
      {
        // we include only contacts that carry current
        Boundary::ConstModelIterator modit(bd->models_begin());
        assert(modit != bd->models_end());

        const DDInterfaceModel* mod = static_cast<const DDInterfaceModel*>(modit->second);
        if (mod->has_current())
        {
          string name(bd->get_name() + ".current");
          // if currents should be plotted, add it also to the plot variables
          if (plot_curr) add_plot_variable(name);
          cnt_set.insert(name);
        }
      }
    }

    // now we declare them
    unsigned int id = static_cast<ID>(ContactCurrent);
    for (set<string>::iterator i(cnt_set.begin()); i != cnt_set.end(); ++i)
    {
      ++id;
      declare_solution_ext(*i, id, SolutionDescriptor::REAL,
          SolutionDescriptor::GLOBAL, units);
    }
  }

  declare_solution(IQE, REAL, GLOBAL, "");
}






void
DriftDiffusion::do_newton(void)
{

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  system.set_options(get_solver_options());
  system.solve();
}








void
DriftDiffusion::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& points)
{

  unsigned int np = points.size();

  ID subdomain = elem->subdomain_id();

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>(0);

  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();
  const libMesh::NumericVector<Number>& oldsolution = system->get_vector("old_sol");

  const unsigned int dim = get_mesh().mesh_dimension();

  const libMesh::DofMap& dof_map = system->get_dof_map();

  unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_e");
  unsigned int ep_var = system->variable_number("fermi_h");
  if (_useparticle == 'e')
    ep_var = en_var;
  else if (_useparticle == 'h')
    en_var = ep_var;

  bool exclude_e = false;
  bool exclude_h = false;
  map<ID, set<string>>::const_iterator it(get_excluded_domains().find(subdomain));
  if (it != get_excluded_domains().end())
  {
    if ((it->second).count("fermi_e"))
      exclude_e = true;
    if ((it->second).count("fermi_h"))
      exclude_h = true;
  }



  libMesh::FEType fe_type = system->variable_type(u_var);
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

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
  libMesh::RealGradient jn(0);
  libMesh::RealGradient jp(0);
  libMesh::RealGradient el_field(0);
  libMesh::RealVectorValue polariz(0);
  double el_pot;

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    double en = 0.0;
    double ep = 0.0;
    double oldu  = 0.0;
    double olden = 0.0;
    double oldep = 0.0;
    double T  = 0.0;
    libMesh::RealGradient e_field(0);
    libMesh::RealGradient grad_en_loc(0);
    libMesh::RealGradient grad_ep_loc(0);
    libMesh::RealGradient grad_T_loc(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices_u[i]);
      en += phi[i][n] * solution(dof_indices_en[i]);
      ep += phi[i][n] * solution(dof_indices_ep[i]);

      oldu += phi[i][n] * oldsolution(dof_indices_u[i]);
      olden += phi[i][n] * oldsolution(dof_indices_en[i]);
      oldep += phi[i][n] * oldsolution(dof_indices_ep[i]);

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

    if (exclude_e) grad_en_loc = 0;
    if (exclude_h) grad_ep_loc = 0;

    el_pot = u;


    sc->set_coordinates(real_pts[n]);

    sc->set_potentials(u, en, ep);
    sc->set_old_potentials(phi0 * oldu, phi0 * olden, phi0 * oldep);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_e(grad_en_loc);
    sc->set_grad_fermi_h(grad_ep_loc);

    sc->calculate_densities();

    double edens = (sc->is_dielectric() ? 0.0 : sc->get_electron_density());
    double hdens = (sc->is_dielectric() ? 0.0 : sc->get_hole_density());
    //double edens = sc->get_electron_density();
    //double hdens = sc->get_hole_density();

    sc->calculate_mobilities();

    sc->compute_thermoelectric_powers();
    double Pn =  sc->get_electron_thermoelectric_power();
    double Pp =  sc->get_hole_thermoelectric_power();

    //sc->calculate_net_recombination_rates();

    double sigma_e = Constants::e * sc->get_electron_conductivity();
    double sigma_h = Constants::e * sc->get_hole_conductivity();

    libMesh::RealGradient dfn = grad_en_loc + Pn * grad_T_loc;
    libMesh::RealGradient dfp = grad_ep_loc + Pp * grad_T_loc;

    libMesh::RealGradient jn_loc = -sigma_e * dfn;
    libMesh::RealGradient jp_loc = -sigma_h * dfp;
    jn += jn_loc;
    jp += jp_loc;

    el_field += e_field;
    polariz += sc->get_total_polarization();


    if (values.count(ElPotential))
      values[ElPotential][n] = u - _reference_potential;

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

    if (values.count(eDensity))
      values[eDensity][n] = edens;

    if (values.count(hDensity))
      values[hDensity][n] = hdens;

    if (values.count(IntrinsicDensity))
      values[IntrinsicDensity][n] = sc->get_intrinsic_density();

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

    bool trapped_electrons = values.count(IonizedElectronTraps);
    bool trapped_holes = values.count(IonizedHoleTraps);
    if (trapped_electrons || trapped_holes)
    {
      sc->calculate_traps();
      if (trapped_electrons)
        values[IonizedElectronTraps][n] = sc->get_ionized_electron_traps();

      if (trapped_holes)
        values[IonizedHoleTraps][n] = sc->get_ionized_hole_traps();
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
      bool get_recomb_e = values.count(eNetRecombination);
      bool get_recomb_h = values.count(hNetRecombination);
      bool get_recomb = get_recomb_e || get_recomb_h;
      double tot_rec_e = 0;
      double tot_rec_h = 0;

      bool need_recomb = get_recomb;
      set<ID>::const_iterator rec_it(_recombination_ids.begin());
      for ( ; rec_it != _recombination_ids.end(); ++rec_it)
      {
        need_recomb |= values.count(*rec_it + eNetRecombination) ||
            values.count(*rec_it + hNetRecombination);
      }

      if (need_recomb)
      {
        // loop over all recombination ids
        rec_it = _recombination_ids.begin();
        for ( ; rec_it != _recombination_ids.end(); ++rec_it)
        {
          bool requested_e = values.count(*rec_it + eNetRecombination);
          bool requested_h = values.count(*rec_it + hNetRecombination);
          pair<double, double> rec;
          if (get_recomb || requested_e || requested_h)
            rec = sc->get_net_recombination_rate(*rec_it);

          if (requested_e)
            values[*rec_it + eNetRecombination][n] = rec.first;
          if (requested_h)
            values[*rec_it + hNetRecombination][n] = rec.second;

          if (get_recomb)
          {
            tot_rec_e += rec.first;
            tot_rec_h += rec.second;
          }
        }

        if (get_recomb_e)
          values[eNetRecombination][n] = tot_rec_e;
        if (get_recomb_h)
          values[hNetRecombination][n] = tot_rec_h;
      }
    }


    {
      bool ept = values.count(ePeltier);
      bool hpt = values.count(hPeltier);

      if (ept || hpt)
        sc->compute_thermoelectric_power_gradient();

      if (ept)
      {
        libMesh::RealGradient PnGrad = sc->get_electron_thermoelectric_power_gradient();
        values[ePeltier][n] =  -T * (PnGrad * jn_loc);
      }

      if (hpt)
      {
        libMesh::RealGradient PpGrad = sc->get_hole_thermoelectric_power_gradient();
        values[hPeltier][n] =  -T * (PpGrad * jp_loc);
      }
    }

    if (values.count(RecombHeat))
    {
      vector<ID> rec_model_ids;
      int n_rec = sc->get_net_recombination_rate_IDs(rec_model_ids);
      double rec_e = 0;
      double rec_h = 0;
      for (int i = 0; i < n_rec; i++)
      {
        pair<double, double> rec(sc->get_net_recombination_rate(rec_model_ids[i]));
        rec_e += rec.first;
        rec_h += rec.second;
      }

      values[RecombHeat][n] = Constants::e * (rec_h * (ep + T * Pp) -
          rec_e * (en + T * Pn));
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

  if (values.count(ElectronBands))
  {
    vector<double> cb;
    sc->get_conduction_bands(cb);
    for (size_t i = 0; i < cb.size(); ++i)
      values[ElectronBands][i] = cb[i] - el_pot;
  }

  if (values.count(HoleBands))
  {
    vector<double> vb;
    sc->get_valence_bands(vb);
    for (size_t i = 0; i < vb.size(); ++i)
      values[HoleBands][i] = vb[i] - el_pot;
  }

}



void
DriftDiffusion::calculate_mean_fermi_levels(void)
{
  ContactData::iterator it = _boundary_eqfermi.begin();
  for ( ; it != _boundary_eqfermi.end(); ++it)
  {
    // zero out current
    (*it).second = 0.0;
  }

  for (it =  _boundary_hqfermi.begin(); it != _boundary_hqfermi.end(); ++it)
  {
    // zero out current
    (*it).second = 0.0;
  }

  ContactData boundary_area;

  SimulationEnvironment& env = get_environment();

  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();
  const libMesh::NumericVector<Number>& solution = system.get_solution_vector();

  const unsigned int dim = mesh.mesh_dimension();

  //const Device& device = *_device;
  //const SimulationEnvironment& environment = get_environment();

  libMesh::DenseVector<Number> X;
  libMesh::DenseSubVector<Number>
                                  Xu(X),
                                  Xn(X),
                                  Xp(X);

  const Options& params = get_my_options();

  // the scaling parameters
  const Scaling& scaling = get_scaling();
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();

  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int en_var = system.variable_number("fermi_e");
  const unsigned int ep_var = system.variable_number("fermi_h");

  FEType fe_type = system.variable_type(en_var);

  libMeshEnums::Order integration_order = params.integration_order;

  // the finite element for boundary integration
  UniquePtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));


  UniquePtr<QBase> qface(QBase::build(
        params.quadrature_type, dim - 1, libMeshEnums::CONSTANT));
  fe_face->attach_quadrature_rule(qface.get());


  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  //
  const vector<vector<libMesh::RealGradient> >&  dphi_face = fe_face->get_dphi();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point_face = fe_face->get_xyz();
  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW_face = fe_face->get_JxW();

  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // we construct a map of all interface models, and remember the
  // associated Boundary pointer.
  map<DDInterfaceModel*, Boundary*> ifmodels;


  BoundaryElementMap::iterator bel(env.boundary_elements_begin());
  const BoundaryElementMap::iterator bend(env.boundary_elements_end());
  //MeshBase::const_element_iterator el =
  //                                this->active_local_elements_begin();
  //const MeshBase::const_element_iterator end_el =
  //                                this->active_local_elements_end();

  // loop over all active elements
  for ( ; bel != bend ; ++bel)
  {
    const Elem* elem = *bel;
    const Elem* top_parent = elem->top_parent();

    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    unsigned int n_dofs     = dof_indices_en.size();
    unsigned int n_dofs_tot = 3 * n_dofs;


    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      Boundary* bd = get_environment().get_boundary(ElementSide(elem, s));
      if (bd == NULL)
        continue;

      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      if (sm != NULL)
      {

        X.resize(n_dofs_tot);    Xu.reposition(0, n_dofs);
        if (_useparticle == 'h')
          Xn.reposition(2 * n_dofs, n_dofs);
        else
          Xn.reposition(n_dofs, n_dofs);
        if (_useparticle == 'e')
          Xp.reposition(n_dofs, n_dofs);
        else
          Xp.reposition(2 * n_dofs, n_dofs);

        dof_map.extract_local_vector(solution, dof_indices, X);

        fe_face->reinit(elem, s);

        int phi_size = phi_face.size();

        // now integrate to include von Neumann and mixed type BCs
        // and polarization
        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {
          // get the solution values at the quadrature point
          Real en = 0.0;
          Real ep = 0.0;
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            en += phi_face[i][qp] * Xn(i);
            ep += phi_face[i][qp] * Xp(i);
          }

          // the jacobian x weight x scaling
          double J = JxW_face[qp];

          _boundary_eqfermi[bd] += J * en * phi0;
          _boundary_hqfermi[bd] += J * ep * phi0;
          boundary_area[bd] += J;
        }

        ifmodels[sm] = bd;

      }
    }
  }

  for (it = _boundary_eqfermi.begin(); it != _boundary_eqfermi.end(); ++it)
  {
    it->second /= boundary_area[it->first];
  }
  for (it = _boundary_hqfermi.begin(); it != _boundary_hqfermi.end(); ++it)
  {
    it->second /= boundary_area[it->first];
  }

  // search for interface models with external current source
  // TODO to be implemented also for holes
  map<const SimulationInterface*, set<Boundary*>> fluxmodels;
  map<DDInterfaceModel*, Boundary*>::iterator ifit(ifmodels.begin());
  for ( ; ifit != ifmodels.end(); ++ifit)
  {
    if ((ifit->first)->get_eflux_simulation() != NULL)
      fluxmodels[(ifit->first)->get_eflux_simulation()].insert(ifit->second);
  }

  map<const SimulationInterface*, double> reffermi_e;
  map<const SimulationInterface*, set<Boundary*>>::iterator sit(fluxmodels.begin());
  for ( ; sit != fluxmodels.end(); ++sit)
  {
    const set<Boundary*>& bdset = sit->second;
    reffermi_e[sit->first] = 0;
    double count = bdset.size();
    for (set<Boundary*>::iterator bdit(bdset.begin()); bdit != bdset.end(); ++bdit)
      reffermi_e[sit->first] += _boundary_eqfermi[*bdit] / count;
  }


  ifit = ifmodels.begin();
  for ( ; ifit != ifmodels.end(); ++ifit)
  {
    double eref = 0.0;
    double href = 0.0;

    if ((ifit->first)->get_eflux_simulation() != NULL)
    {
      map<const SimulationInterface*, double>::iterator simit(
          reffermi_e.find((ifit->first)->get_eflux_simulation()));
      if (simit != reffermi_e.end())
        eref = simit->second;
    }

    (ifit->first)->set_reference_fermi_potentials(eref, href);
  }
}


void
DriftDiffusion::calculate_currents_rstf_global(void)
{

  if (_rstf == NULL)
    prepare_rstf();

  // now we have certainly the RSTFs prepared

  // put all RSTFs into a map
  map<const Boundary*, libMesh::NumericVector<double>*> rstf;

  {
    ContactData::iterator it = _boundary_currents.begin();
    for ( ; it != _boundary_currents.end(); ++it)
    {
      // zero out current
      (*it).second = 0.0;

      rstf[it->first] = _rstf->get_testfunction(it->first);
    }
  }


  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>(0);

  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  //const Device& device = *(_device);
  SimulationEnvironment& env = get_environment();

  const libMesh::DofMap& dof_map = system->get_dof_map();
  const libMesh::DofMap& dof_map_rstf = _rstf->get_dof_map();

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

  libMesh::FEType fe_type = system->variable_type(u_var);

  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
  libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(
        get_my_options().quadrature_type, dim, get_my_options().integration_order));
  fe->attach_quadrature_rule(qrule.get());


  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe->get_JxW();
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  // physical coordinates of the quadrature points
  //const vector<Point>& q_point = fe->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();


  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  vector<unsigned int> dof_indices_rstf;

  MeshBase::const_element_iterator el(this->active_local_elements_begin());
  const MeshBase::const_element_iterator end_el(this->active_local_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    //const Elem* top_parent = elem->top_parent();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    dof_map_rstf.dof_indices(elem, dof_indices_rstf, 0);

    DDBulkModel* sc =
        get_bulk_model<DDBulkModel>(elem);


    assert(sc != NULL);


    // in a dielectric we have no current
    if (sc->is_dielectric())
      continue;


    fe->reinit(elem);

    sc->reinit(elem);

    //Get the temperature given the element
    vector<double> T_nodes =  sc->get_temperature_at_nodes();


    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      unsigned int n_dofs = dof_indices_u.size();
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      libMesh::RealGradient dEfn(0);
      libMesh::RealGradient dEfp(0);
      libMesh::RealGradient e_field(0);
      libMesh::RealGradient dT(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        ep += phi[i][qp] * solution(dof_indices_ep[i]);

        dEfn += dphi[i][qp] * solution(dof_indices_en[i]);
        dEfp += dphi[i][qp] * solution(dof_indices_ep[i]);

        dT += dphi[i][qp] * T_nodes[i];

        e_field -= dphi[i][qp] * solution(dof_indices_u[i]);
      }

      // prepare for calculating local properties
      //sc->set_coordinates(elem->centroid());  ????? 2012-08-31
      sc->set_coordinates(q_point[qp]);


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
      sc->set_electric_field(phi0 * e_field);
      sc->set_grad_fermi_e(phi0 * dEfn);
      sc->set_grad_fermi_h(phi0 * dEfp);

      sc->calculate_densities();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      //Get the thermoelectric power
      sc->compute_thermoelectric_powers();
      double Pn =  sc->get_electron_thermoelectric_power() / phi0;
      double Pp =  sc->get_hole_thermoelectric_power() / phi0;

      // we put the minus here for convenience
      double sigma_e = -Constants::e * sc->get_electron_conductivity();
      double sigma_h = -Constants::e * sc->get_hole_conductivity();

      double Rn = sc->get_net_electron_recombination_rate();
      double Rp = sc->get_net_hole_recombination_rate();
      double net_rate = JxW[qp] * Constants::e * (Rn - Rp);

      // loop over all recombination models and add the ones that model
      // tunneling for the current contact
      // this is needed, because the current calculated afterwards does not
      // include direct tunneling.
      DriftDiffusionProperties::RecombinationModelIterator recit(
          sc->recombination_models_begin());
      DriftDiffusionProperties::RecombinationModelIterator recend(
          sc->recombination_models_end());
      for ( ; recit != recend; ++recit)
      {
        const Boundary* bd = (*recit)->get_tunneling_contact();
        if (bd != NULL)
        {
          double elrate, hlrate;
          (*recit)->get_net_recombination_rates(elrate, hlrate);
          _boundary_currents[bd] -= JxW[qp] * Constants::e * (elrate - hlrate);
        }
      }

      libMesh::RealGradient je(JxW[qp] * phi0 * (sigma_e * (dEfn + Pn * dT)));
      libMesh::RealGradient jh(JxW[qp] * phi0 * (sigma_h * (dEfp + Pp * dT)));

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        // do this for each contact
        map<const Boundary*, libMesh::NumericVector<double>*>::iterator rstf_it(rstf.begin());
        for ( ; rstf_it != rstf.end(); ++rstf_it)
        {
          const Boundary* bd = rstf_it->first;
          const libMesh::NumericVector<double>& sol = *rstf_it->second;

          _boundary_currents[bd] += ((je + jh) * dphi[n][qp] +
              net_rate * phi[n][qp]) * sol(dof_indices_rstf[n]);

        }
      }

    } // end loop over quadrature points
  } // end loop over elements

}



void
DriftDiffusion::calculate_currents_rstf_compact(void)
{

  // reset currents
  {
    ContactData::iterator it =
      _boundary_currents.begin();
    for ( ; it != _boundary_currents.end(); ++it)
      (*it).second = 0.0;
  }

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>(0);

  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  //const Device& device = *(_device);
  SimulationEnvironment& env = get_environment();

  const libMesh::DofMap& dof_map = system->get_dof_map();

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

  libMesh::FEType fe_type = system->variable_type(u_var);

  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
  libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(
        get_my_options().quadrature_type, dim, get_my_options().integration_order));
  fe->attach_quadrature_rule(qrule.get());


  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe->get_JxW();
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  // physical coordinates of the quadrature points
  //const vector<Point>& q_point = fe->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();


  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;


  BoundaryElementMap::iterator el(env.boundary_elements_begin());
  BoundaryElementMap::iterator end_el(env.boundary_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    //const Elem* top_parent = elem->top_parent();

    const Boundary* boundary = el.get_boundary();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DDBulkModel* sc =
        get_bulk_model<DDBulkModel>(elem);

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
      Real u  = 0.0;
      Real en = 0.0;
      Real ep = 0.0;
      libMesh::RealGradient dEfn(0);
      libMesh::RealGradient dEfp(0);
      libMesh::RealGradient e_field(0);
      libMesh::RealGradient dT(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        ep += phi[i][qp] * solution(dof_indices_ep[i]);

        dEfn += dphi[i][qp] * solution(dof_indices_en[i]);
        dEfp += dphi[i][qp] * solution(dof_indices_ep[i]);

        dT += dphi[i][qp] * T_nodes[i];

        e_field -= dphi[i][qp] * solution(dof_indices_u[i]);
      }

      // prepare for calculating local properties
      //sc->set_coordinates(elem->centroid());  ????? 2012-08-31
      sc->set_coordinates(q_point[qp]);


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
      sc->set_electric_field(phi0 * e_field);
      sc->set_grad_fermi_e(phi0 * dEfn);
      sc->set_grad_fermi_h(phi0 * dEfp);

      sc->calculate_densities();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      //Get the thermoelectric power
      sc->compute_thermoelectric_powers();
      double Pn =  sc->get_electron_thermoelectric_power() / phi0;
      double Pp =  sc->get_hole_thermoelectric_power() / phi0;

      // we put the minus here for convenience
      double sigma_e = -Constants::e * sc->get_electron_conductivity();
      double sigma_h = -Constants::e * sc->get_hole_conductivity();

      double Rn = sc->get_net_electron_recombination_rate();
      double Rp = sc->get_net_hole_recombination_rate();
      double net_rate = JxW[qp] * Constants::e * (Rn - Rp);

      libMesh::RealGradient je(JxW[qp] * phi0 * (sigma_e * (dEfn + Pn * dT)));
      libMesh::RealGradient jh(JxW[qp] * phi0 * (sigma_h * (dEfp + Pp * dT)));

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        _boundary_currents[boundary] += ((je + jh) * dphi[n][qp] -
            net_rate * phi[n][qp]) * weight[n];

    } // end loop over quadrature points
  } // end loop over elements

}



void
DriftDiffusion::calculate_field_emission(void)
{

  const SimulationEnvironment& env = get_environment();

  ContactData fe_currents;

  const set<PhysicalModel*>& pm = get_interface_models();
  set<PhysicalModel*>::const_iterator it(pm.begin());
  const set<PhysicalModel*>::const_iterator end(pm.end());
  for ( ; it != end; ++it)
  {
    DDInterfaceModel* mod = static_cast<DDInterfaceModel*>(*it);

    // register the contact if it is a real contact (with current)
    if (mod->has_field_emission())
    {
      const Boundary* bnd = get_environment().get_boundary(mod->get_name());
      if (bnd != NULL)
      {
        fe_currents[bnd] = 0.0;
      }
    }
  }

  // if there is nothing to do we return immediately
  if (fe_currents.empty()) return;

  // otherwise we have to calculate it

  string file = get_option("field_emission_file", "");

  bool write_emission = file.empty() ? false : true;

  ofstream of;
  if (write_emission)
  {
    file = get_output_directory() + "/" + file;
    of.open(file.c_str());
    if (!of.good())
      throw RuntimeException("Cannot open file " + file + " for writing.");

    of << "% Field emission currents calculated by TiberCAD\n"
       << "% units: SI\n"
       << "%\n"
       << "% Columns:\n"
       << "% pos_x  pos_y  pos_z  mom_x  mom_y  mom_z  mass  charge  current\n";
  }


  TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = get_mesh();

  const libMesh::DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = get_scaling().get_potential_scaling();


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");

  libMesh::FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = libMeshEnums::FIRST;

  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
        get_my_options().quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe_face->get_JxW();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<libMesh::RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;


  // Note: the following is not very general at the moment

  // we only loop over the boundary sides


  MeshBase::const_element_iterator el =
                                  this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    //ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);

    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

    assert(sc != NULL);

    // for now, we calculate field emission only in dielectrics
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

        DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

        // check if we should do something
        if (!sm->has_field_emission())
          break;

        sc->reinit(elem);

        fe_face->reinit(elem, s);

        int phi_size = phi.size();

        double current = 0.0;

        FowlerNordheim* em = sm->get_field_emission_model();
        sm->reinit(elem, s);

        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {
          sm->set_face_normal(face_normals[qp]);

          // get the solution value at the quadrature point
          libMesh::RealGradient e_field(0.0);
          for (unsigned int i = 0; i < phi_size; i++)
            e_field += dphi[i][qp] * solution(dof_indices_u[i]);

          sc->set_electric_field(phi0 * e_field);

          double F = phi0 * e_field * face_normals[qp];
          double curr = em->get_emission_current(F);
          current += JxW[qp] * curr;

        } // end loop over quadrature points

        if (write_emission)
        {
          double v = em->get_velocity() / 100.0;
          v /= (Constants::c);
          v /= -sqrt(1.0 - v * v); // need inward direction
          libMesh::RealVectorValue(q_point[0] *
              get_scaling().get_calc_mesh_units()).write_unformatted(of, false);
          libMesh::RealVectorValue(face_normals[0] * v).write_unformatted(of, false);
          of << "9.11e-31 -1.6e-19 " << current << "\n";
        }

        fe_currents[boundary] += current;
      }
    } // end loop over elem sides
  } // end loop over elements


  Messages m;
  m.info("Field emission currents:");
  m.indent();
  for (ContactData::iterator it(fe_currents.begin()); it != fe_currents.end(); ++it)
  {
    ostringstream os;
    os << (*it).first->get_name() << ": " << (*it).second << "" << endl;
    m.info(os.str());
  }
}




void
DriftDiffusion::calculate_currents_surfint(void)
{

  // reset currents
  ContactData::iterator it =
    _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second = 0.0;

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const libMesh::NumericVector<Number>& solution = get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const libMesh::DofMap& dof_map = system->get_dof_map();

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

  libMesh::FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = get_my_options().integration_order;

  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
        get_my_options().quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe_face->get_JxW();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<libMesh::RealGradient> >& dphi = fe_face->get_dphi();

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

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

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
            libMesh::RealGradient e_field(0);
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
          libMesh::RealGradient e_field(0);
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
          sc->set_grad_fermi_e(phi0 * libMesh::RealGradient(dEfn, 0, 0));
          sc->set_grad_fermi_h(phi0 * libMesh::RealGradient(dEfp, 0, 0));

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
DriftDiffusion::calculate_surface_recombination(void)
{

  // we only do something if we are on processor 0
  //if (get_communicator().rank() != 0)
  //  return;

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const libMesh::NumericVector<Number>& solution = get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const libMesh::DofMap& dof_map = system->get_dof_map();

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

  libMesh::FEType fe_type = system->variable_type(u_var);

  // the finite element for boundary integration
  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type));
  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order = get_my_options().integration_order;

  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
        get_my_options().quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());


  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe_face->get_JxW();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<libMesh::RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_ep;

  // the total recombination current
  double current = 0.0;

  //BoundaryElementMap::iterator el(env.boundary_elements_begin());
  //BoundaryElementMap::iterator end_el(env.boundary_elements_end());
  MeshBase::const_element_iterator el =
                                  this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_en, en_var);
    dof_map.dof_indices(elem, dof_indices_ep, ep_var);

    //DDBulkModel* sc =
    //  dynamic_cast<DDBulkModel*>(get_physical_model(subdomain));

    //assert(sc != NULL);

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      // the check about the current is a bit primitive, but ok for now
      if ((sm != NULL)  && !sm->has_current())
      {
        sm->reinit(elem, s);

        fe_face->reinit(elem, s);

        int phi_size = phi.size();

        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {
          // get the solution value at the quadrature point
          Real u  = 0.0;
          Real en = 0.0;
          Real ep = 0.0;
          Real dEfn = 0.0;
          Real dEfp = 0.0;
          libMesh::RealGradient e_field(0);
          for (unsigned int i = 0; i < phi_size; i++)
          {
            u  += phi[i][qp] * solution(dof_indices_u[i]);
            en += phi[i][qp] * solution(dof_indices_en[i]);
            ep += phi[i][qp] * solution(dof_indices_ep[i]);

            double tmp = dphi[i][qp] * face_normals[qp];
            dEfn += tmp * solution(dof_indices_en[i]);
            dEfp += tmp * solution(dof_indices_ep[i]);

            e_field += dphi[i][qp] * solution(dof_indices_u[i]);
          }

          // prepare for calculating local properties
          sm->set_coordinates(q_point[qp]);
          sm->set_potentials(phi0 * u, phi0 * en, phi0 * ep);
          sm->set_electric_field(-phi0 * e_field);
          sm->set_grad_fermi_e(phi0 * dEfn);
          sm->set_grad_fermi_h(phi0 * dEfp);
          sm->set_face_normal(face_normals[qp]);
          sm->compute();

          const vector<double>& coeff_a = sm->get_a();
          const vector<double>& coeff_g = sm->get_g();

          //double tmp = (coeff_g[0] - coeff_a[0] * u * phi0);
          //if (sm->get_type(0) != DDInterfaceModel::DIRICHLET)
          //  value_u += J * tmp / (x0 * C0);
          //else
          //  value_u = coeff_g[0] / phi0;

          double value_n = 0;
          if (sm->get_type(1) != DDInterfaceModel::DIRICHLET)
            value_n = (coeff_g[1] - coeff_a[1] * en);


          double value_p = 0;
          if (sm->get_type(2) != DDInterfaceModel::DIRICHLET)
            value_p = (coeff_g[2] - coeff_a[2] * ep);


          current += JxW[qp] * 0.5 * (value_n + value_p);

        } // end loop over quadrature points
      }
    } // end loop over elem sides
  } // end loop over elements

  // accumulate from all
  this->get_communicator().sum(current);

  ostringstream rec;
  rec << "Surface recombination current = " << current * Constants::e << "\n";
  Messages::info(rec.str());
}



void
DriftDiffusion::build_local_scaling(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  const libMesh::NumericVector<Number>& solution = get_solution_vector();
  libMesh::NumericVector<Number>& loc_scaling = system.get_vector("scaling");
  loc_scaling.zero();

  if (!_do_local_scaling)
  {
    // we do not use local scaling, therefore set it to 1
    loc_scaling.add(1.0);
    loc_scaling.close();
    return;
  }

  // aliases for nicer code
  const MeshBase& mesh = get_mesh();

  const libMesh::DofMap& dof_map = system.get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const Options& params = get_my_options();

  // the scaling parameters to scale back the result
  const Scaling& scaling = get_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double x0 = scaling.get_length_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;
  // density scaling for electrons
  double C0_e = C0;
  // density scaling for holes
  double C0_h = C0;
  // scaling for recombination rates
  double R0_e = C0_e / scaling.get_time_scaling();
  double R0_h = C0_h / scaling.get_time_scaling();

  C0_e = C0_h = 1;

  const unsigned int u_var = system.variable_number("potential");
  // fermi potential variable numbers are defined within element loop
  // since (in future) different variables can ben set in different regions

  // we have to detect the dirichlet DOFs
  set<unsigned int> dirichlet_dofs;

  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  // dof indices vectors for fermi potentials are defined within element loop
  // since (in future) different variables can ben set in different regions

  libMesh::FEType fe_type = system.variable_type(u_var);
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type, true));
  libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(
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
  const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();



  //
  // for boundary elements
  //

  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));
  libMeshEnums::Order integration_order = params.integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
        params.quadrature_type, dim - 1, integration_order));
  fe_face->attach_quadrature_rule(qface.get());

  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  //
  const vector<vector<libMesh::RealGradient> >&  dphi_face = fe_face->get_dphi();
  //
  // physical coordinates of the quadrature points
  const vector<Point>& q_point_face = fe_face->get_xyz();
  //
  const vector<Point>& face_normals = fe_face->get_normals();
  //
  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW_face = fe_face->get_JxW();

  libMesh::DenseVector<Number> local_scaling;
  libMesh::DenseSubVector<Number> scaleu(local_scaling);
  // scale vectors for fermi potentials are defined within element loop
  // since (in future) different variables can ben set in different regions

  MeshBase::const_element_iterator it =
    this->active_local_elements_begin();
  const MeshBase::const_element_iterator end =
    this->active_local_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);
    assert(sc != NULL);

    sc->reinit(elem);
    fe->reinit(elem);

    //Get variables for fermi potentials
    map<unsigned int, string> q_var;
    for (auto&& cp : sc->get_carrier_properties())
      q_var.insert( make_pair(system.variable_number(cp.first), cp.first) );

    //Get dof indices
    map<unsigned int, vector<unsigned int>> dof_indices_q;

    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    for (auto var : q_var)
    {
      vector<unsigned int> dof_indices_q_tmp;
      dof_map.dof_indices(elem, dof_indices_q_tmp, var.first);
      dof_indices_q.insert( make_pair(var.first, dof_indices_q_tmp) );
    }

    unsigned int n_dofs = dof_indices_u.size();

    local_scaling.resize(dof_indices.size());
    scaleu.reposition(0, n_dofs);
    //Define scale vector for fermi potentials
    map<unsigned int, libMesh::DenseSubVector<Number>> scaleq;
    unsigned int n_var = 0;
    for (auto var : q_var)
    {
      n_var++;
      libMesh::DenseSubVector<Number> scale_tmp(local_scaling);
      scale_tmp.reposition(n_var*n_dofs, n_dofs);
      scaleq.insert( make_pair(var.first, scale_tmp) );
    }

    assert(elem->n_nodes() == dof_indices_u.size());

    libMesh::RealGradient field(0.0);
    for (unsigned int i = 0; i < dof_indices_u.size(); i++)
      field += dphi[i][0] * solution(dof_indices_u[i]);
    field *= -phi0;

    
    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {
      // get the solution values at the quadrature point
      Real u  = 0.0;
      map<unsigned int, Real> q;
      for (auto var : q_var)
        q.insert( make_pair(var.first, 0.0) );

      libMesh::RealGradient e_field(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        for (auto var : q_var)
          q[var.first] += phi[i][qp] * solution(dof_indices_q[var.first][i]);
        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }


      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_el_potential(phi0 * u);
      for (auto var : q_var)
      {
        sc->set_fermi_potential(var.second, phi0 * q[var.first]);
        sc->set_grad_fermi(var.second, libMesh::RealGradient(0.0));
      }
      sc->set_electric_field(phi0 / x0 * e_field);

      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();

      const libMesh::RealTensor& permittivity = sc->get_relative_permittivity();

      double l2_eps = JxW[qp] * l2;

      map<unsigned int, double> sigma, ddens_dphi;

      for (auto var : q_var)
      {
        sigma.insert( make_pair(var.first, JxW[qp] * sc->get_q_conductivity(var.second) / mu0 / C0) );
        ddens_dphi.insert( make_pair(var.first, sc->get_q_density_derivative(var.second)) );
      }

      /*
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
      */

      double drhovec[2];
      sc->get_charge_density_derivatives(drhovec);
      double drho = -JxW[qp] * (drhovec[0] + drhovec[1]) * phi0 / C0;
      if (sc->is_dielectric())
        drho = 0.0;
      if (params.local_neutrality)
        drho = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        double phi_i_x_phi_j = JxW[qp] * phi[i][qp] * phi[i][qp];

        for (auto& scale : scaleq)
          scale.second(i) += sigma[scale.first] * (dphi[i][qp] * dphi[i][qp]);
                             //- dRn[1] * phi_i_x_phi_j;

        scaleu(i) +=
            //phi_i_x_phi_j;
	  l2_eps * (dphi[i][qp] * (permittivity * dphi[i][qp])) -
            drho * phi[i][qp] * phi[i][qp];
      }

    } // end loop over quadrature points



    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      if (sm != NULL)
      {
        for (unsigned int i = 0; i < n_dofs; i++)
        {
          if (elem->is_node_on_side(i, s))
          {
            if (sm->get_type(0) == DDInterfaceModel::DIRICHLET)
              dirichlet_dofs.insert(dof_indices_u[i]);
              //scaleu(i) *= _penalty_value;

            n_var = 0;
            for (auto& dofs : dof_indices_q)
            {
              n_var++;
              if (sm->get_type(n_var) == DDInterfaceModel::DIRICHLET)
              dirichlet_dofs.insert(dofs.second[i]);
            }

            //if (sm->get_type(1) == DDInterfaceModel::DIRICHLET)
              //dirichlet_dofs.insert(dof_indices_en[i]);
                //scalen(i) *= _penalty_value;
          }
        }

      }
    }



    /*
    // loop over sides
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      bool true_boundary = get_environment().is_outer_boundary(ElementSide(elem, s));

      if ((sm != NULL) || true_boundary)
      {
        fe_face->reinit(elem, s);

        int phi_size = phi_face.size();


        // now integrate to include von Neumann and mixed type BCs
        // and polarization
        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {

          double epsilon = sc->get_relative_permittivity();
          double l2_eps = l2 * epsilon;


          // get the solution values at the quadrature point
          Real u  = 0.0;
          Real en = 0.0;
          Real ep = 0.0;
          libMesh::RealGradient e_field(0);
          libMesh::RealGradient grad_en(0);
          libMesh::RealGradient grad_ep(0);
          libMesh::RealGradient grad_T(0);
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            u  += phi_face[i][qp] * solution(dof_indices_u[i]);
            en += phi_face[i][qp] * solution(dof_indices_en[i]);
            ep += phi_face[i][qp] * solution(dof_indices_ep[i]);
            e_field -= dphi_face[i][qp] * solution(dof_indices_u[i]);
            grad_en += dphi_face[i][qp] * solution(dof_indices_en[i]);
            grad_ep += dphi_face[i][qp] * solution(dof_indices_ep[i]);
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

          / *
          sc->compute_thermoelectric_powers();
          double Pn =  sc->get_electron_thermoelectric_power() / phi0;
          double Pp =  sc->get_hole_thermoelectric_power() / phi0;

          double jn = 0.0;
          double jp = 0.0;
          if (coupling & ECURRENT)
            jn = (sigma_e * grad_en + Pn * grad_T) * face_normals[qp] / x0;
          if (coupling & HCURRENT)
            jp = -(sigma_h * grad_ep + Pp * grad_T) * face_normals[qp] / x0;
          * /

          // the jacobian x weight x scaling
          double J = JxW_face[qp];


          // contribution to the jacobian
          if (sm != NULL)
          {
            sm->compute();

            double scale_u = J * phi0 / x0 / C0;
            if (sm->get_type(0) == DDInterfaceModel::DIRICHLET)
              scale_u *= _penalty_value;

            double scale_n = J* phi0 / (x0 * R0_e);
            if (sm->get_type(1) == DDInterfaceModel::DIRICHLET)
              scale_n *= _penalty_value;

            double scale_p = J * phi0 / (x0 * R0_h);
            if (sm->get_type(2) == DDInterfaceModel::DIRICHLET)
              scale_p *= _penalty_value;


            const vector<double>& deriv_u = sm->get_jacobian_row(0);
            const vector<double>& deriv_en = sm->get_jacobian_row(1);
            const vector<double>& deriv_ep = sm->get_jacobian_row(2);


            for (unsigned int i = 0; i < n_dofs; i++)
            {

              Real phi_i_x_phi_j =
                phi_face[i][qp] * phi_face[i][qp];

              local_scaling_[elem->get_node(i)][2] -= scale_u * deriv_u[0] * phi_i_x_phi_j;


              ////
              // NOTE:
              //   the signs are inverted here because outflow means recombination
              //   and the normal point outwards, giving a positiv current for
              //   outflow

              local_scaling_[elem->get_node(i)][0] -= scale_n * deriv_en[1] * phi_i_x_phi_j;

              local_scaling_[elem->get_node(i)][1] += scale_p * deriv_ep[2] * phi_i_x_phi_j;
            }
          }
        }
      }
    }
    */


    loc_scaling.add_vector(local_scaling, dof_indices);
  } // end loop over elements
  loc_scaling.close();

  set<unsigned int>::iterator dirit(dirichlet_dofs.begin());
  const set<unsigned int>::iterator dirend(dirichlet_dofs.end());
  for ( ; dirit != dirend; ++dirit)
  {
    double val = loc_scaling.el(*dirit) * _penalty_value;
    loc_scaling.set(*dirit, val);
  }

  loc_scaling.close();


  /*
  {
    map<const Node*, vector<double> >::iterator it(local_scaling_.begin());
    map<const Node*, vector<double> >::iterator end(local_scaling_.end());
    for ( ; it != end; ++it)
    {
      (it->second)[0] = sqrt((it->second)[0]);
      (it->second)[1] = sqrt((it->second)[1]);
      (it->second)[2] = sqrt((it->second)[2]);
    }
  }
  */
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

    if (id == IQE)
    {
      values[id] = vector<double>(1, _iqe);
    }
    else
    {
      // for now it can only be currents
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
}



void
DriftDiffusion::calculate_currents(void)
{
  if (get_my_options().current_calculation == RSTF)
    calculate_currents_rstf_global();
  else if (get_my_options().current_calculation == RSTF_COMPACT)
    calculate_currents_rstf_compact();
  else
    calculate_currents_surfint();


  // sum up contributions from all processes
  ContactData::iterator it = _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
  {
    this->get_communicator().sum((*it).second);
  }
}





double
DriftDiffusion::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_potential_scaling();
}





void
DriftDiffusion::assemble_system(const libMesh::NumericVector<Number>& x,
    libMesh::NumericVector<Number>* residual,
    libMesh::SparseMatrix<Number>* jacobian,
    libMesh::NonlinearImplicitSystem&)
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
      break;
  }

}






template <int coupling>
void
DriftDiffusion::do_assembly(const libMesh::NumericVector<Number>& x,
    libMesh::NumericVector<Number>* residual,
    libMesh::SparseMatrix<Number>* jacobian)
{
  START_LOG(get_name() + ": Matrix assembly", "");

  //build_local_scaling();

  // references for nicer code
  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const unsigned int dim = mesh.mesh_dimension();

  //const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  const Options& params = get_my_options();


  libMesh::NumericVector<Number>& oldx = system.get_vector("old_sol");
  libMesh::NumericVector<Number>& loc_scaling = system.get_vector("scaling");
  //SparseMatrix<double>& sysmat = system.get_matrix("Preconditioner");
  //if (residual != NULL)
  //if (jacobian != NULL)
  //  sysmat.zero();

  //
  // some scaling stuff...
  //
  // NOTE: the mesh and all parameters were not explicitly scaled, so
  //       we have to treat scaling by explicit division/multiplication
  //
  // the scaling parameters
  const Scaling& scaling = get_scaling();
  // the scaling parameter for the poisson eq.
  // The factor 1e-2 comes from the fact, that we are calculating in cm!
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double C0 = (_do_local_scaling) ? 1.0 : scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  // x 1e4 because we calculate in cm, but P comes in C/m^2
  const double P0 = (Constants::e * x0 * C0) * 1e4;
  // density scaling for electrons
  double C0_e = C0;
  // density scaling for holes
  double C0_h = C0;


  // scaling for recombination rates
  double R0 = C0 / scaling.get_time_scaling();

  const libMesh::DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  // fermi potential variable numbers are defined within element loop
  // since (in future) different variables can ben set in different regions

  libMesh::FEType fe_type = system.variable_type(u_var);

  libMeshEnums::Order integration_order = params.integration_order;

  // the finite element
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type, true));
  libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(
        params.quadrature_type, dim, integration_order));
  fe->attach_quadrature_rule(qrule.get());

  // the finite element for boundary integration
  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));

  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;

  libMesh::UniquePtr<libMesh::QBase> qface(libMesh::QBase::build(
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
  const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();



  // references to boundary-specific data that will be used to
  // assemble the system.
  // Data will be given for each quadrature point.
  //
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  //
  const vector<vector<libMesh::RealGradient> >&  dphi_face = fe_face->get_dphi();
  //
  // physical coordinates of the quadrature points
  const vector<libMesh::Point>& q_point_face = fe_face->get_xyz();
  //
  const vector<libMesh::Point>& face_normals = fe_face->get_normals();
  //
  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW_face = fe_face->get_JxW();


  // the system matrix (will hold also element jacobian contribution)
  libMesh::DenseMatrix<Number> Ke;
  // the system rhs (will hold also element rhs contribution)
  libMesh::DenseVector<Number> Fe;
  // the local solution
  libMesh::DenseVector<Number> X;
  // the local old step
  libMesh::DenseVector<Number> oldX;
  // the local scaling
  DenseVector<Number> local_scaling;


  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  // dof indices vectors for fermi potentials are defined within element loop
  // since (in future) different variables can ben set in different regions

  // zero out residual and jacobian !! IMPORTANT !!
  if (residual != NULL)
    residual->zero();
  if (jacobian != NULL)
    jacobian->zero();


  MeshBase::const_element_iterator el =
                                  this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  this->active_local_elements_end();

  // loop over all active elements
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    ID subdomain = elem->subdomain_id();

    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);

    assert(sc != NULL);
    sc->reinit(elem);

    //Get variables for fermi potentials
    map<unsigned int, string> q_var;
    map<unsigned int, unsigned int> q_var_bc;
    q_var.insert( make_pair(u_var, "potential") );
    for (auto&& cp : sc->get_carrier_properties())  //get carriers
      q_var.insert( make_pair(system.variable_number(cp.first), cp.first) ); 
                            //insert variable number   and      variable name (the same as carrier name)
    for (auto var : q_var)
      q_var_bc.insert( make_pair(var.first, q_var_bc.size()) );

   //Get dof indices
    map<unsigned int, vector<unsigned int>> dof_indices_var;

    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    for (auto var : q_var)
    {
      vector<unsigned int> dof_indices_var_tmp;
      dof_map.dof_indices(elem, dof_indices_var_tmp, var.first);
      dof_indices_var.insert( make_pair(var.first, dof_indices_var_tmp) ); 
                            //insert variable number and dof indices for the element
    }

    unsigned int n_dofs     = dof_indices_u.size();
    unsigned int n_dofs_tot = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs_tot, n_dofs_tot);
    Fe.resize(n_dofs_tot);
    X.resize(n_dofs_tot);
    oldX.resize(n_dofs_tot);
    local_scaling.resize(n_dofs_tot);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);
    dof_map.extract_local_vector(oldx, dof_indices, oldX);
    dof_map.extract_local_vector(loc_scaling, dof_indices, local_scaling);

    //define submatrices
    map<unsigned int, map<unsigned int, libMesh::DenseSubMatrix<Number>>> Kvv;
    map<unsigned int, libMesh::DenseSubVector<Number>> Fv,
                                                       Xv,
                                                       oldXv,
                                                       scalev;

    // Reposition the submatrices according to this scheme:
    //
    //        | Kuu Kun Kup Ku. |        | Fu |
    //   Ke = | Knu Knn Knp ... |;  Fe = | Fn |
    //        | Kpu Kpn Kpp ... |        | Fp |
    //        | K.u ... ... ... |        | ...|
    //

    unsigned int n_var = 0;
    for (auto var : q_var)
    {
      libMesh::DenseSubVector<Number> F_tmp(Fe), 
                                      X_tmp(X), 
                                      oldX_tmp(oldX), 
                                      scale_tmp(local_scaling);

      F_tmp.reposition(n_var*n_dofs, n_dofs);
      X_tmp.reposition(n_var*n_dofs, n_dofs);
      oldX_tmp.reposition(n_var*n_dofs, n_dofs);
      scale_tmp.reposition(n_var*n_dofs, n_dofs);

      Fv.insert( make_pair(var.first, F_tmp) );
      Xv.insert( make_pair(var.first, X_tmp) );
      oldXv.insert( make_pair(var.first, oldX_tmp) );
      scalev.insert( make_pair(var.first, scale_tmp) );

      map<unsigned int, libMesh::DenseSubMatrix<Number>> Kv_tmp; //Kvv columns
      unsigned int n_varj = 0;
      for (auto varj : q_var)
      {
        libMesh::DenseSubMatrix<Number> K_tmp(Ke);
        K_tmp.reposition(n_var*n_dofs, n_varj*n_dofs, n_dofs, n_dofs);
        Kv_tmp.insert( make_pair(varj.first, K_tmp) );
        n_varj++;
      }
      Kvv.insert( make_pair(var.first, Kv_tmp) );
      n_var++;
    }

    // Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_at_nodes();

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      map<unsigned int, Real> u;
      map<unsigned int, Real> oldu;
      map<unsigned int, libMesh::RealGradient> grad_u;

      for (auto var : q_var)
      {
        u.insert( make_pair(var.first, 0.0) );
        oldu.insert( make_pair(var.first, 0.0) );
        grad_u.insert( make_pair(var.first, libMesh::RealGradient(0)) );
      }

      // get the solution values at the quadrature point
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (auto var : q_var)
        {
          u[var.first] += phi[i][qp] * Xv[var.first](i);
          oldu[var.first] += phi[i][qp] * oldXv[var.first](i);
          grad_u[var.first] += dphi[i][qp] * Xv[var.first](i);
        }
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      double grad_fac = phi0 / x0;
      for (auto var : q_var)
      {
        if (var.first == u_var)
        {
          sc->set_el_potential(phi0 * u[var.first]);
          sc->set_old_el_potential(phi0 * oldu[var.first]);
          sc->set_electric_field(-grad_fac * grad_u[var.first]);
        }
        else
        {
          sc->set_fermi_potential(var.second, phi0 * u[var.first]);
          sc->set_old_fermi_potential(var.second, phi0 * oldu[var.first]);
          sc->set_grad_fermi(var.second, grad_fac * grad_u[var.first]);
        }
      }

      // calculate all local properties
      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      map<unsigned int, long double> dens, R;
      map<unsigned int, double> mu, sigma;

      for (auto var : q_var)
      {
        if (var.first != u_var)
        {
          dens.insert( make_pair(var.first, sc->get_q_density(var.second)) );
          R.insert( make_pair(var.first, sc->get_net_q_recombination_rate(var.second)) );
          mu.insert( make_pair(var.first, sc->get_q_mobility(var.second)) );
          sigma.insert( make_pair(var.first, sc->get_q_conductivity(var.second) / (mu0 * C0)) );
        }
      }
      //double Nd = sc->get_ionized_donor_density();
      //double Na = sc->get_ionized_acceptor_density();

      /*
      // Get the thermoelectric power
      sc->compute_thermoelectric_powers();
      double eTEpower =  sc->get_electron_thermoelectric_power() / phi0;
      double hTEpower =  sc->get_hole_thermoelectric_power() / phi0;
      */

      const libMesh::RealTensor& permittivity = sc->get_relative_permittivity();

      // the jacobian x weight x scaling
      double J = JxW[qp];

      /*
      double sigma_e_x_Pe_x_J = J * sigma_e * eTEpower;
      double sigma_h_x_Ph_x_J = J * sigma_h * hTEpower;
      */

      // TEST
      //double dn_dphi = sc->get_electron_density_derivative();
      //double art_diff = 0.5 * x0 * elem->hmax() * mue * dn_dphi * grad_en.size() / (mu0 * C0_e);
      //cerr << "art. diffusivity = " << art_diff << " (" << sigma_e << ")" << endl;
      //sigma_e += art_diff;

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
          double laplace =
	    J * (dphi[i][qp] * dphi[j][qp]);

	  double laplace_u =
            J * (dphi[i][qp] * (permittivity * dphi[j][qp]));

          for (auto var : q_var)
          {
            if (var.first == u_var)
            {
              if (coupling & POISSON)
                Kvv[var.first][var.first](i,j) += l2 * laplace_u / scalev[var.first](i);
              else
                Kvv[var.first][var.first](i,j) += 1.0;
            }
            else
            {
              if (coupling & CURRENTS)
                Kvv[var.first][var.first](i,j) += sigma[var.first] * laplace / scalev[var.first](i);
              else
                Kvv[var.first][var.first](i,j) += 1.0;
            }


          }
        }


      }

      //
      // for jacobian compute the other contributions
      //
      if (jacobian != NULL)
      {
        map<unsigned int, long double> ddens_dphi;
        map<unsigned int, map<unsigned int, long double>> dR_ddens,
                                                          dR_dEf,
                                                          dR;

        for (auto vari : q_var)
        {
          if (vari.first != u_var)
          {
            long double ddens_dphi_tmp = sc->get_q_density_derivative(vari.second);
            ddens_dphi.insert( make_pair(vari.first, ddens_dphi_tmp) );
            dR[vari.first].insert( make_pair(u_var, 0.0) );

            for (auto varj : q_var)
            {
              if (varj.first != u_var)
              {
                long double dR_ddens_tmp = sc->get_net_q_recombination_rate_derivatives(vari.second, varj.second)[0];
                long double dR_dEf_tmp = sc->get_net_q_recombination_rate_derivatives(vari.second, varj.second)[1];
                long double dR_tmp = (dR_dEf_tmp - dR_ddens_tmp * ddens_dphi_tmp) * phi0 / R0;

                if (R[vari.first] == 0.0)
                  dR_ddens_tmp = dR_dEf_tmp = dR_tmp = 0.0;

                dR_ddens[vari.first].insert( make_pair(varj.first, dR_ddens_tmp) );
                dR_dEf[vari.first].insert( make_pair(varj.first, dR_dEf_tmp) );
                dR[vari.first].insert( make_pair(varj.first, dR_tmp) );
                dR[vari.first][u_var] -= dR_ddens_tmp * ddens_dphi_tmp * phi0 / R0;
              }

            }
          }
        }

        map<unsigned int, double> drho;
        drho.insert( make_pair(u_var, 0.0) );
        //drho[1] = (dn_dphi - dNd_dphi) * phi0 / C0;
        //drho[2] = -(dp_dphi - dNa_dphi) * phi0 / C0;
        for (auto var : q_var)
        {
          if (var.first != u_var)
          {
            char ctype = sc->get_carrier_properties(var.second)->get_carrier_type();
            drho[var.first] = (ctype == 'e') ? ddens_dphi[var.first] * phi0 / C0 : -ddens_dphi[var.first] * phi0 / C0;
            drho[u_var] -= drho[var.first];

            if (params.local_neutrality)
              drho[u_var] = drho[var.first] = 0.0;
            if (sc->is_dielectric())
              drho[u_var] = drho[var.first] = 0.0;
          }
        }


        // d(sigma_n)/du * element-jacobian
        // sigma_n = mu_n * n means the conductivity of electrons
        // the factor phi_0 comes from the derivative with respect to the potentia
        map<unsigned int, double> dsigma, dmu_u;
        map<unsigned int, libMesh::RealGradient> dmu_grad_u;
        map<unsigned int, libMesh::RealGradient> dmu_grad_f;

        for (auto var : q_var)
        {
          if (var.first != u_var)
          {
            double dsigma_tmp = J * phi0 / (mu0 * C0_e) * mu[var.first] * ddens_dphi[var.first];
            dsigma.insert( make_pair(var.first, dsigma_tmp) );

            double dmu_u_tmp = sc->get_q_mobility_derivative_potential(var.second);
            dmu_u_tmp *= J * phi0 / (mu0 * C0) * dens[var.first];
            dmu_u.insert( make_pair(var.first, dmu_u_tmp) );

            libMesh::RealGradient dmu_grad_u_tmp(0);
            sc->get_q_mobility_derivative_grad_potential(var.second, dmu_grad_u_tmp);
            dmu_grad_u_tmp *= J * phi0 / (mu0 * C0) * dens[var.first] / x0;
            dmu_grad_u.insert( make_pair(var.first, dmu_grad_u_tmp) );

            if (dim > 1)
            {
              libMesh::RealGradient dmu_grad_f_tmp(0);
              sc->get_q_mobility_derivative_grad_fermi(var.second, dmu_grad_f_tmp);
              dmu_grad_f_tmp *= J * phi0 / (mu0 * C0) * dens[var.first] / x0;
              dmu_grad_f.insert( make_pair(var.first, dmu_grad_f_tmp) );
            }
          }
        }

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          map<unsigned int, double> lap,
                                    dsigma_x_lap;
          for (auto var : q_var)
          {
            lap.insert( make_pair(var.first, (dphi[i][qp] * grad_u[var.first]) / scalev[var.first](i)) );
            dsigma_x_lap.insert ( make_pair(var.first, dsigma[var.first]*lap[var.first]) );
          }

          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            //

            if (_options.exact_newton)
            {

              for (auto var : q_var)
              {
                if ((var.first != u_var) && (coupling & CURRENTS))
                {
                  double dsigma_x_phi = dsigma_x_lap[var.first] * phi[j][qp];
                  double dmu_u_x_phi = dmu_u[var.first] * phi[j][qp];
                  double dmu_grad_u_x_dphi = dmu_grad_u[var.first] * dphi[j][qp];
                  double dmu_grad_f_x_dphi = dmu_grad_f[var.first] * dphi[j][qp];

                  if (coupling & POISSON)
                    Kvv[var.first][u_var](i,j) += dsigma_x_phi + (dmu_u_x_phi + dmu_grad_u_x_dphi) * lap[var.first];

                  Kvv[var.first][var.first](i,j) += (dmu_grad_f_x_dphi - dmu_u_x_phi) * lap[var.first] - dsigma_x_phi;
                }
              }
            }


            /*
            // contribution of the Seebeck effect ->residual_derivative
            double dsigma_e_x_phi_x_Pe = dsigma_e * phi[j][qp] * eTEpower;
            double dsigma_h_x_phi_x_Ph = dsigma_h * phi[j][qp] * hTEpower;

            for (unsigned int k = 0; k < n_dofs; k++)
            {
              double laplace = dphi[i][qp] * dphi[k][qp];

              if (coupling & ECURRENT)
              {
                double elem_contrib =
                  dsigma_e_x_phi_x_Pe * laplace * T_nodes[k] / scalen(i);

		Knn(i,j) -= elem_contrib;

                if (coupling & POISSON)
		  Knu(i,j) += elem_contrib;

              }

              if (coupling & HCURRENT)
              {
                double elem_contrib =
                  dsigma_h_x_phi_x_Ph * laplace * T_nodes[k] / scalep(i);

		Kpp(i,j) -= elem_contrib;

		if (coupling & POISSON)
		  Kpu(i,j) += elem_contrib;

		}
            }
            */



            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (coupling & POISSON)
            {
              Kvv[u_var][u_var](i,j) -= drho[u_var] * phi_i_x_phi_j / scalev[u_var](i);

              for (auto var : q_var)
              {
                if ( (var.first != u_var) && (coupling & CURRENTS) )
                  Kvv[u_var][var.first](i,j) -= drho[var.first] * phi_i_x_phi_j / scalev[u_var](i);
              }
            }

            if (coupling & CURRENTS)
            {
              for (auto vari : q_var)
              {
                if (vari.first != u_var)
                {
                  if ((coupling & POISSON) && _options.exact_newton)
                    Kvv[vari.first][u_var](i,j) -= dR[vari.first][u_var] * phi_i_x_phi_j / scalev[vari.first](i);

                  for (auto varj : q_var)
                    Kvv[vari.first][varj.first](i,j) -= dR[vari.first][varj.first] * phi_i_x_phi_j / scalev[vari.first](i);
                }
              }
            }

          }
        }

      } //end jacobian


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        long double J_x_rho;
        J_x_rho = J * sc->get_charge_density() / C0;


        if (params.local_neutrality)
          J_x_rho = 0.0;
        if (sc->is_dielectric())
          J_x_rho = 0.0;

        long double J_x_P0 = J / P0;

        // net recombination rate
        map<unsigned int, long double> J_x_R;
        for (auto var : q_var)
        {
          if (var.first != u_var)
            J_x_R.insert( make_pair(var.first, J*R[var.first]/R0) );
        }

        /*
        libMesh::RealVectorValue P(sc->get_total_polarization());
        P *= J_x_P0;
        */

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (auto var : q_var)
          {
            if (var.first == u_var)
            {
              if (coupling & POISSON)
                Fv[var.first](i) -= J_x_rho * phi[i][qp] / scalev[var.first](i);
              else
                Fv[var.first](i) -= Xv[var.first](i);
            }
            else
            {
              long double net_recomb = J_x_R[var.first] * phi[i][qp] / scalev[var.first](i);

              if (coupling & CURRENTS)
                Fv[var.first](i) -= net_recomb;
              else
                Fv[var.first](i) -= Xv[var.first](i);
            }

          }
        }

        /*
        // include Seebeck contribution -> Residual
	for (unsigned int i = 0; i < n_dofs; i++)
	{
	  for (unsigned int k = 0; k < n_dofs; k++)
	  {
	    Real laplace = dphi[i][qp] * dphi[k][qp];

	    if (coupling & ECURRENT)
	      Fn(i) += sigma_e_x_Pe_x_J * laplace * T_nodes[k] / scalen(i);

	    if (coupling & HCURRENT)
	      Fp(i) += sigma_h_x_Ph_x_J * laplace * T_nodes[k] / scalep(i);
	  }
	}
        */

      }  //end residual

    } // end loop over quadrature points



    map<unsigned int, vector<double> > dirichlet_jac;
    map<unsigned int, double> dirichlet_res;

    // now loop over the element sides to find boundary elements
    // and to include von Neumann and mixed type boundary conditions
    //
    // NOTE 1:
    // we dont apply BC for nabla(Ef) but for the particle
    // flux mu * n * nabla(Ef)
    //
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {

      Material* mat = get_material(elem);
      DDInterfaceModel* sm = get_interface_model<DDInterfaceModel>(elem, s);

      bool true_boundary = environment.is_outer_boundary(ElementSide(elem, s));

      if (sm != NULL)
      {
        fe_face->reinit(elem, s);

        sm->reinit(elem, s);

        int phi_size = phi_face.size();

        // now integrate to include von Neumann and mixed type BCs
        // and polarization
        for (unsigned int qp = 0; qp < qface->n_points(); qp++)
        {

          map<unsigned int, Real> u;
          map<unsigned int, Real> oldu;
          map<unsigned int, libMesh::RealGradient> grad_u;

          for (auto var : q_var)
          {
            u.insert( make_pair(var.first, 0) );
            grad_u.insert( make_pair(var.first, libMesh::RealGradient(0)) );
          }

          // get the solution values at the quadrature point
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (auto var : q_var)
            {
              u[var.first] += phi_face[i][qp] * Xv[var.first](i);
              grad_u[var.first] += dphi_face[i][qp] * Xv[var.first](i);
            }
          }

          // prepare for calculating local properties
          sc->set_coordinates(q_point_face[qp]);

          double grad_fac = phi0 / x0;
          for (auto var : q_var)
          {
            if (var.first == u_var)
            {
              sc->set_el_potential(phi0 * u[var.first]);
              sc->set_electric_field(-grad_fac * grad_u[var.first]);
            }
            else
            {
              sc->set_fermi_potential(var.second, phi0 * u[var.first]);
              sc->set_grad_fermi(var.second, grad_fac * grad_u[var.first]);
            }
          }

          sc->calculate_densities();
          sc->calculate_mobilities();

          // we put the phi0 here for convenience
          //double sigma_e = phi0 * sc->get_electron_conductivity();
          //double sigma_h = phi0 * sc->get_hole_conductivity();

          //sc->compute_thermoelectric_powers();
          //double Pn =  sc->get_electron_thermoelectric_power() / phi0;
          //double Pp =  sc->get_hole_thermoelectric_power() / phi0;

          //double jn = 0.0;
          //double jp = 0.0;
          //if (coupling & ECURRENT)
          //  jn = (sigma_e * grad_en + Pn * grad_T) * face_normals[qp] / x0;
          //if (coupling & HCURRENT)
          //  jp = -(sigma_h * grad_ep + Pp * grad_T) * face_normals[qp] / x0;

          sm->set_face_normal(face_normals[qp]);
          sm->set_coordinates(q_point_face[qp]);

          for (auto var : q_var)
          {
            if (var.first == u_var)
            {
              sm->set_el_potential(phi0 * u[var.first]);
              sm->set_electric_field(-grad_fac * grad_u[var.first]);
            }
            else
            {
              sm->set_fermi_potential(var.second, phi0 * u[var.first]);
              sm->set_grad_fermi(var.second, grad_fac * grad_u[var.first]);
            }
          }

          sm->compute();

          // the jacobian x weight x scaling
          double J = JxW_face[qp];


          // contribution to the jacobian
          if (jacobian != NULL)
          {
            // for Dirichlet DOFs we do not add anything
            vector<double> scale_v(q_var_bc.size());
            vector<vector<double>> deriv_v(q_var_bc.size());
            for (auto var : q_var_bc)
            {
              scale_v[var.second] = (var.first==u_var) ?  J * phi0 / x0 / C0 : J * phi0 / (x0 * R0);
              if (sm->get_type(var.second) == DDInterfaceModel::DIRICHLET)
                scale_v[var.second] = 0;

              vector<double> deriv_tmp = sm->get_jacobian_row(var.second);
              deriv_v[var.second] = deriv_tmp;
            }


            for (unsigned int i = 0; i < n_dofs; i++)
            {
              vector<double> fac(q_var_bc.size());
              for (auto var : q_var_bc)
              {
                fac[var.second] = scale_v[var.second] / scalev[var.first](i);

                if (elem->is_node_on_side(i, s))
                {
                  if (sm->get_type(var.second) == DDInterfaceModel::DIRICHLET)
                    dirichlet_jac[var.second*n_dofs + i] = deriv_v[var.second];
                }
              }

              for (unsigned int j = 0; j < n_dofs; j++)
              {
                Real phi_i_x_phi_j = phi_face[i][qp] * phi_face[j][qp];

                for (auto vari : q_var_bc)
                {
                  if (((vari.first == u_var) && (coupling & POISSON)) || ((vari.first != u_var) && (coupling & CURRENTS)))
                  {
                    for (auto varj : q_var_bc)
                    {
                      if (((varj.first == u_var) && (coupling & POISSON)) || ((varj.first != u_var) && (coupling & CURRENTS)))
                        Kvv[vari.first][varj.first](i,j) -= fac[vari.second] * deriv_v[vari.second][varj.second] * phi_i_x_phi_j;
                    }
                  }
                }
                ////
                // NOTE:
                //   the signs are inverted here because outflow means recombination
                //   and the normal point outwards, giving a positive current for
                //   outflow
              }
            }
          } // end jacobian

          // contribution to -Fe_i
          if (residual != NULL)
          {

            vector<double> value_v(q_var_bc.size());

            const vector<double>& coeff_a = sm->get_a();
            const vector<double>& coeff_g = sm->get_g();

            for (auto var : q_var_bc)
            {
              value_v[var.second] = (var.first == u_var) ? J / (x0 * C0) : J / (x0 * R0);
              if (sm->get_type(0) != DDInterfaceModel::DIRICHLET)
                value_v[var.second] *= (coeff_g[var.second] - coeff_a[var.second] * u[var.first] * phi0);
              else
                value_v[var.second] = coeff_g[var.second] / phi0;
            }

            /*
            if (0) //(true_boundary)
            {
              // If we are on an outer boundary, we include
              // the polarization

              if (params.default_boundary_condition == ZEROFIELD)
              {

                libMesh::RealVectorValue P(0.0);
                P = sc->get_total_polarization();
                double Pn = (P * face_normals[qp]) / P0;
                double value_u = -J * Pn;

                if (coupling & POISSON)
                {
                  for (unsigned int i = 0; i < n_dofs; i++)
                  {
                    Fu(i) -= value_u * phi_face[i][qp] / scaleu(i);
                  }
                }
              }
            }
            */

            for (unsigned int i = 0; i < n_dofs; i++)
            {
              for (auto var : q_var_bc)
              {
                double scale_v = 1.0 / scalev[var.first](i);

                if (elem->is_node_on_side(i, s))
                {
                  if (sm->get_type(var.second) == DDInterfaceModel::DIRICHLET)
                  {
                    dirichlet_res[var.second*n_dofs + i] = value_v[var.second];
                    scale_v = 0;
                  }
                }

                if (((var.first == u_var) && (coupling & POISSON)) || ((var.first != u_var) && (coupling & CURRENTS)))
                  Fv[var.first](i) -= value_v[var.second] * phi_face[i][qp] * scale_v;

              } //end q_var_bc loop
            }  // end dof loop
          }  // end residual
        } // end loop over qp
      } // end if (sm != NULL)
      else if ((true_boundary) && (residual != NULL))
      {
        // If we are on an outer boundary, we have to include
        // the polarization
        //
        // NOTE:
        // we only include the polarization when no explicit
        // boundary is defined

        if (params.default_boundary_condition == ZEROFIELD)
        {
          fe_face->reinit(elem, s);

          for (unsigned int qp = 0; qp < qface->n_points(); qp++)
          {
            /*
            // the jacobian x weight x scaling
            double J = JxW_face[qp];

            libMesh::RealVectorValue P(0.0);
            P = sc->get_total_polarization();
            double Pn = (P * face_normals[qp]) / P0;
            double value_u = -J * Pn;

            if (coupling & POISSON)
            {
              for (unsigned int i = 0; i < n_dofs; i++)
              {
                Fu(i) -= value_u * phi_face[i][qp] / scaleu(i);
              }
            }
            */

          }
        }
      }  // end true boundary && residual
    } // end loop over element side


    // impose Dirichlet BCs
    if (jacobian != NULL)
    {
      map<unsigned int, vector<double> >::iterator it(dirichlet_jac.begin());
      const map<unsigned int, vector<double> >::iterator end(dirichlet_jac.end());
      for ( ; it != end; ++it)
      {
        unsigned int i = it->first;
        for (unsigned int j = 0; j < n_dofs_tot; ++j)
          Ke(i,j) = 0.0;

        unsigned int j = i % n_dofs;
        for (auto var : q_var_bc)
          Ke(i,var.second * n_dofs + j) = -(it->second)[var.second];
      }
    }

    if (residual != NULL)
    {
      map<unsigned int, double>::iterator it(dirichlet_res.begin());
      const map<unsigned int, double>::iterator end(dirichlet_res.end());
      for ( ; it != end; ++it)
      {
        unsigned int i = it->first;
        for (unsigned int j = 0; j < n_dofs_tot; ++j)
          Ke(i,j) = 0.0;
        Ke(i,i) = 1;

        Fe(i) = -(it->second);
      }
    }

    // check if it is a dielectric
    if (sc->is_dielectric())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (auto vari : q_var)
        {
          for (auto varj : q_var)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
              Kvv[vari.first][varj.first](i,j) = ((vari.first == u_var) && (varj.first == u_var)) ? Kvv[vari.first][varj.first](i,j) : 0.0;
          }

          if (vari.first != u_var)
          {
            if (is_dielectric_boundary_node(elem->get_node(i)))
              Kvv[vari.first][vari.first](i, i) = 0.0;
            else
              Kvv[vari.first][vari.first](i, i) = 1.0;
            // we simply set the electrochemical potential to zero
            Fv[vari.first](i) = 0.0;
          }
        }
      }
    }


    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    // NOTE: this changes dof_indices that's why the application of
    //       Dirichlet type BCs needs special care
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    system.exclude_dofs(Ke, dof_indices, elem);


    //perf_log.start_event("add");

    if (residual != NULL)
    {
      for (unsigned int i = 0; i < n_dofs_tot; i++)
        for (unsigned int j = 0; j < n_dofs_tot; j++)
          Fe(i) += Ke(i,j) * x(dof_indices[j]);

      if (coupling & ELECTRONS)
      {
        //cerr << Ke << endl;
        //TiberMath::svd(Ke, Fe);
        //elem->centroid().write_unformatted(cerr, false);
        //cerr << " " << Fe(0) <<  "  " << Fe(Fe.size()-1) << endl;
      }

      system.exclude_dofs(Fe, dof_indices, elem);

      residual->add_vector(Fe, dof_indices);
      //sysmat.add_matrix(Ke, dof_indices);
    }
    else
    {

      jacobian->add_matrix(Ke, dof_indices);
      /*
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        unsigned int i2 = i + n_dofs;
        unsigned int i3 = i2 + n_dofs;
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          if (j != i)
          {
            Ke(i, j) = Ke(i, j + n_dofs) = Ke(i, j + 2*n_dofs) = 0.0;
            Ke(i2, j) = Ke(i2, j + n_dofs) = Ke(i2, j + 2*n_dofs) = 0.0;
            Ke(i3, j) = Ke(i3, j + n_dofs) = Ke(i3, j + 2*n_dofs) = 0.0;
          }
        }
      }
      sysmat.add_matrix(Ke, dof_indices);
      */
      //cerr << Ke << endl;
      //TiberMath::svd(Ke, Fe);
      //elem->centroid().write_unformatted(cerr, false);
      //cerr << " " << Fe(0) <<  "  " << Fe(Fe.size()-1) << endl;
    }


  } // end loop over elements

  if (jacobian != NULL)
  {
    jacobian->close();

    /*
    DenseMatrix<Number> Pe(3, 3);
    Pe.zero();

    unsigned int sysnum = system.get_libmesh_system()->number();
    MeshBase::const_node_iterator nit(this->active_nodes_begin());
    const MeshBase::const_node_iterator nend(this->active_nodes_end());
    for ( ; nit != nend; ++nit)
    {
      if ((*nit)->has_dofs(sysnum))
      {
        unsigned int dof_1 = (*nit)->dof_number(sysnum, 0, 0);
        unsigned int dof_2 = (*nit)->dof_number(sysnum, 1, 0);
        unsigned int dof_3 = (*nit)->dof_number(sysnum, 2, 0);

        double a11 = (*jacobian)(dof_1, dof_1);
        double a12 = (*jacobian)(dof_1, dof_2);
        double a13 = (*jacobian)(dof_1, dof_3);
        double a21 = (*jacobian)(dof_2, dof_1);
        double a22 = (*jacobian)(dof_2, dof_2);
        double a23 = (*jacobian)(dof_2, dof_3);
        double a31 = (*jacobian)(dof_3, dof_1);
        double a32 = (*jacobian)(dof_3, dof_2);
        double a33 = (*jacobian)(dof_3, dof_3);

        Pe(0, 0) = a22*a33 - a32*a23;
        Pe(0, 1) = a13*a32 - a33*a12;
        Pe(0, 2) = a12*a23 - a22*a13;
        Pe(1, 0) = a23*a31 - a33*a21;
        Pe(1, 1) = a11*a33 - a31*a13;
        Pe(1, 2) = a13*a21 - a23*a11;
        Pe(2, 0) = a21*a32 - a31*a22;
        Pe(2, 1) = a12*a31 - a32*a11;
        Pe(2, 2) = a11*a22 - a21*a12;
        double det = a11 * Pe(0, 0) + a21 * Pe(0, 1) + a31 * Pe(0, 2);
        Pe.scale(1.0 / det);

        sysmat.add_matrix(Pe, {dof_1, dof_2, dof_3});



      }
    }
    sysmat.close();
    */
    /*
    if (coupling & ELECTRONS)
    {
      ostringstream os;
      os << "_" << __private_counter << ".m";
      ostringstream ms;
      ms << "writing " << "J" << os.str() << "\n";
      Messages::info(ms.str());
      jacobian->print_matlab("J" + os.str());
      sysmat.print_matlab("precond" + os.str());
    }
    if (coupling & ELECTRONS) __private_counter++;
    //if (__private_counter == 2) exit(0);
    */
  }
  else
  {
    residual->close();

    //sysmat.close();
    //sysmat.print_matlab("sysmat.m");
    /*
    if (coupling & ELECTRONS)
    {
      //if (__private_counter > 0)
      //  oldx -= x;
      ostringstream os;
      os << "_" << __private_counter;
      write_nodal_vector("residual" + os.str(), *residual);
      ostringstream ms;
      ms << "writing " << "residual" << os.str() << " (norm = " << residual->l2_norm() << ")\n";
      Messages::info(ms.str());
      //residual->print_matlab("F" + os.str() + ".m");
      //write_nodal_vector("x" + os.str(), oldx);
      //NumericVector<Number>& weight = system.get_vector("scaling");
      //write_nodal_vector("scaling" + os.str(), weight);
      __private_counter++;
      //oldx = x;
    }
    */
  }


  STOP_LOG(get_name() + ": Matrix assembly", "");
}





void
DriftDiffusion::do_load_data(istream& is)
{
  SimulationInterface::do_load_data(is);

  compute_scaling(get_my_options().scaling_type);

  build_local_scaling();
}




void
DriftDiffusion::write_nodal_vector(const string& filename, const libMesh::NumericVector<double>& vec)
{

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const libMesh::DofMap& dof_map = system->get_dof_map();

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
    this->active_local_elements_begin();
  const MeshBase::const_element_iterator end =
    this->active_local_elements_end();

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

  vector<string> format;
  get_output_format(format);
  DataOutput data_output(get_mesh(), format[0]);
  data_output.set_output_directory(get_output_directory());
  vector<string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  data_output.write_nodal_data(filename, results, names);

}

