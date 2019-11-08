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
#include "libmesh/node.h"
#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/perf_log.h"
//#include "libMeshDefs.h"

#include "DataOutput.h"
#include "Messages.h"

// C++ includes
#include <fstream>
#include <numeric>

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
    _iqe(0.0),
    _reference_potential(0.0),
    _rstf(NULL)
{

  // first we detect the carriers declared, because we need
  // this in an early stage

  auto physit(get_options().submodels_begin("Physics"));
  ModelOptions& physopts = physit->second;

  auto itc (physopts.submodels_begin("carrier"));
  auto end_itc (physopts.submodels_end("carrier"));

  // names can appear several times
  set<string> names;
  // check if carriers with the same name have also the same charge
  map<string, set<double>> charges;
   // check if carriers with the same name have also the same spin
  map<string, set<double>> spin;

  for ( ; itc != end_itc; ++itc)
  {
    string name = (itc->second).get_option("name", "");
    names.insert(name);
    charges.insert( make_pair(name, set<double>()) );
    charges[name].insert( (itc->second).get_option("charge", 1.0) );
    spin.insert( make_pair(name, set<double>()) );
    spin[name].insert( (itc->second).get_option("spin", 0.5) );
  }

  _carriers.resize(0);
  _carriers.reserve(names.size());
  for (auto&& name : names)
  {
    if (charges[name].size()>1)
      throw InitFailedException("Carrier '" + name + "' has multiple definitions with different charges");

    if (spin[name].size()>1)
      throw InitFailedException("Carrier '" + name + "' has multiple definitions with different spins");

    _carriers.push_back(name);
  }


  // detect all recombination models in order to define specific plot variables
  for (auto&& block : set<string>({"recombination", "trap", "generation"}))
  {

  auto itr (physopts.submodels_begin(block));
  auto end_itr (physopts.submodels_end(block));

  for ( ; itr != end_itr; ++itr)
  {
    string recname = (itr->second).get_option("name", "");
    string plotname = (itr->second).get_option("plot_name", "");
    vector<string> rec_carriers;

    // get the carriers associated to each recombination model
    (itr->second).get_option("carriers", rec_carriers);

    if (plotname != "")
    {
      _rec_models.insert(make_pair(plotname, set<unsigned int>()));

      for (auto carrier : rec_carriers)
      {
        for (unsigned int i = 0; i < _carriers.size(); i++)
        {
          if (carrier == _carriers[i])
            _rec_models[plotname].insert(i);
        }
      }
    }
  } //end rec models
  }

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

  if (model == nullptr)
    throw ModelErrorException(
        "DriftDiffusion: No such physical model: " + modelname);

  // set the ordered list of carriers
  model->set_known_carriers(_carriers);

  return(model);
}



PhysicalModel*
DriftDiffusion::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{

  DDInterfaceModel* model = DDInterfaceModel::create(boundary, options);

  // set the ordered list of carriers
  if (model != nullptr)
    model->set_known_carriers(_carriers);

  return(model);
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

    for (auto&& cp : sc->get_carrier_properties())
    {
      double mu = sc->get_q_mobility(cp.first);
      mu0 = (mu0 > mu) ? mu0 : mu;
    }

    // I don't know what is better...
    double C = fabs(sc->get_material()->get_net_doping_density());
    //double C = fabs(sc->get_ionized_donor_density() -
    //    sc->get_ionized_acceptor_density());
    C0 = (C0 > C) ? C0 : C;

    //sc->get_bulk_equilibrium_densities(densities);
    //double ni = sc->get_intrinsic_density();
    //ni0 = (ni0 > ni) ? ni0 : ni;

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
      //ni0 = 1e12;
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

  libMesh::NumericVector<Number>& solution = system.get_local_solution_vector();

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

  solution.close();
  system.update();
}




void
DriftDiffusion::set_hole_fermi_level(double Ef_p)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = system.get_local_solution_vector();

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

  solution.close();
  system.update();
}




void
DriftDiffusion::set_electric_potential(double pot)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = system.get_local_solution_vector();

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

  solution.close();
  system.update();
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
DriftDiffusion::calculate_conserved_carriers(void)
{
  if (!_conservation.empty())
  {
    for (auto&& cons : _conservation)
    {
      cons.second.volume = 0.0;
      if (cons.second.initial_density > 0)
        cons.second.conserved_number = 0.0;
    }



    MeshBase::const_element_iterator el =
        this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
        this->active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
      DDBulkModel* sc = get_bulk_model<DDBulkModel>(*el);

      // Get variables for fermi potentials
      set<ID> q_var;

      // we will loop only over carriers that are present in this element
      for (auto&& cp : sc->get_carrier_properties())
        q_var.insert(cp.first);


      for (auto&& cons : _conservation)
      {
        bool is_present = false;
        for (auto&& var : q_var)
        {
          if (find(cons.second.carrier_vars.begin(),
              cons.second.carrier_vars.end(), var) !=
                  cons.second.carrier_vars.end())
          {
            is_present = true;
            break;
          }
        }

        if (is_present)
        {
          cons.second.volume += (*el)->volume();
        }
      }
    }

    double scale = 100 * this->get_mesh_units();
    switch (this->get_mesh().mesh_dimension())
    {
      case 3:
        scale *= 100 * this->get_mesh_units();

      case 2:
        scale *= 100 * this->get_mesh_units();

      default:
        break;
    }

    vector<double> vols;
    vols.reserve(_conservation.size());

    for (auto&& cons : _conservation)
    {
      vols.push_back(cons.second.volume);
    }

    // in parallel run need to sum up all values
    this->get_communicator().sum(vols);

    int idx = 0;
    for (auto&& cons : _conservation)
    {
      cons.second.volume = vols[idx];
      if (cons.second.initial_density > 0)
      {
        cons.second.conserved_number = vols[idx] * scale * cons.second.initial_density;
      }
      else
      {
        cons.second.initial_density = cons.second.conserved_number / (vols[idx] * scale);
      }
      idx++;
    }

  }

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
          this->get_environment().get_boundary(it->first)->models_begin()->second);

      if (bd != NULL)
      {
        double voltage = bd->get_simulation_voltage();

        _voltages[it->first] = voltage;

        if (voltage != 0.0)
          equilibrium = false;
      }
    }
  }


  calculate_conserved_carriers();



  if (!equilibrium_done())
  {
    solve_equilibrium();

    //build_local_scaling();

    // if we would repeat the equilibrium simulation, we can stop now
    if (equilibrium)
      return;
  }

  // set the old solution
  //EquationSystems& es = get_equation_systems();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();
  system.get_solution_vector().close();
  system.get_vector("old_sol") = system.get_solution_vector();


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
    do_newton();
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
  //build_local_scaling();

  get_my_options().coupling = coupling;

  // calculate the currents to print them on screen
  calculate_currents();
  calculate_iqe();
  calculate_surface_recombination();

  if (!_conservation.empty())
  {
    Messages m;
    m.newline();
    m.info("Number conservation:");
    m.indent();

    for (auto&& cons : _conservation)
    {

      ostringstream os;
      os << "carriers: ";
      for (auto&& var : cons.second.carrier_vars)
        os << system.variable_name(var) << ", ";
      m.info(os.str());

      const libMesh::DofMap& dof_map = system.get_dof_map();
      vector<dof_id_type> scalars;
      dof_map.SCALAR_dof_indices(scalars, cons.first);
      os.str("");
      const double phi0 = get_scaling().get_potential_scaling();
      os << "mean electrochemical potential: " <<
          phi0 * (system.get_solution_vector())(scalars[0]);
      m.info(os.str());

      os.str("");
      os << "total number of carriers      : " <<
          cons.second.conserved_number;
      m.info(os.str());
    }
    m.newline();
  }

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
    Boundary* bnd = this->get_environment().get_boundary(it->first);

    if (bnd == nullptr)
      continue;

    DDInterfaceModel* ifmod =
        static_cast<DDInterfaceModel*>(bnd->models_begin()->second);

    // we print only contacts with a current
    if (!ifmod->has_current())
      continue;

    ostringstream os;
    os << setprecision(6);
    os << it->first;
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
    os << it->second * bnd->get_area_factor();

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

}


void
DriftDiffusion::calculate_weights(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = system.get_solution_vector();
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
    map<ID, string> var_q;
    for (auto&& cp : sc->get_carrier_properties())
      var_q.insert( make_pair(system.variable_number(_carriers[cp.first]), _carriers[cp.first]) );

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
        sc->set_fermi_potential(var.first, phi0 * q[var.first]);
        sc->set_old_fermi_potential(var.first, phi0 * q[var.first]);
      }

      sc->calculate_densities();

      weight.set(dofu, 1);

      for (auto var : var_q)
      {
        double dens = min(1.0, log10(1 + sc->get_q_density(var.first)) / 18);// < 1e7 ? 0 : 1;
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
      map<ID, vector<double>> eb;
      map<ID, vector<double>> hb;
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

    //declare_solution(ElectronBands, NTUPLE, CELL, "eV", num_eb);
    //declare_solution(HoleBands, NTUPLE, CELL, "eV", num_hb);
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
    const Boundary* bd = this->get_environment().get_boundary(it->first);
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
        sim_voltages[bd->get_name()] = cnt->get_simulation_voltage();
        cnt->set_simulation_voltage(0.0);
      }
    }
  }

  // make a rough guess
  guess_equilibrium();

  //build_local_scaling();


  bool limitstep = get_option("limit_step", false);
  get_options().set_option("limit_step", false);
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
  //build_local_scaling();
  get_options().set_option("limit_step", limitstep);

  // set the contact voltages back to the desired values
  it = sim_voltages.begin();
  for ( ; it != end; ++it)
  {
    const Boundary* bd = this->get_environment().get_boundary(it->first);
    vector<ID> ids;
    bd->get_region_ids(ids);
    for (unsigned int i = 0; i < ids.size(); i++)
    {
      // it must exist
      MaterialBoundary* mb = get_environment().get_device().get_boundary_object(ids[i]);
      ElectricalContact* cnt =
          dynamic_cast<ElectricalContact*>(mb->get_model(get_id()));

      if (cnt != NULL)
        cnt->set_simulation_voltage(sim_voltages[bd->get_name()]);
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

      const libMesh::NumericVector<Number>& solution = system.get_solution_vector();
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

  ID rec_id = INVALID_ID;
  set<ID> nonrad_id;

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

    string radname = (opts->second).get_option("radiative_recombination", "");
    rec_id = get_solution_id(radname);

    vector<string> nradname;
    (opts->second).get_option("nonradiative_recombination", nradname);
    for (auto&& nonrad : nradname)
      nonrad_id.insert(get_solution_id(nonrad));

  }

  if ((rec_id == INVALID_ID) ||
      nonrad_id.count(INVALID_ID))
  {
    if (plot_solution(IQE))
    {
      Messages::warning(string("Cannot calculate IQE: ") +
          "you need to specify explicitly radiative_recombination " +
          "nonradiative_recombination.");
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

  // to accumulate the total nonradiative rec.
  double Rnrad = 0.0;

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
    vector<double>& data = datamap[rec_id];

    for (auto&& nrad : nonrad_id)
      datamap[nrad] = vector<double>(qrule->n_points());

    DriftDiffusion::get_solution_secure(elem, datamap, qrule->get_points());

    double iqe_el = 0;

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {
      iqe_el += JxW[qp] * data[qp];
      for (auto&& nrad : nonrad_id)
        Rnrad += JxW[qp] * datamap[nrad][qp];
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

  vector<double> recs = {_iqe, Rnrad};
  this->get_communicator().sum(recs);
  _iqe = recs[0];
  Rnrad = recs[1];

  double Rtot = _iqe + Rnrad;

  ostringstream rec;
  rec << "Rrad = " << _iqe * Constants::e << 
    " Rnrad = " << Rnrad * Constants::e << "\n";
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
      const Boundary* bnd = this->get_environment().get_boundary(it->first);
      if (bnd == nullptr)
        continue;

      DDInterfaceModel* ifmod =
          static_cast<DDInterfaceModel*>(bnd->models_begin()->second);

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

  libMesh::NumericVector<Number>& solution_u = poisson.get_local_solution_vector();
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

  // for conserved particles, we put a chemical potential based on the
  // prescribed density and the volume
//  for (auto&& scalar : _conservation)
//  {
//    unsigned int dof = scalar.first;
//    vector<dof_id_type> dof_ids;
//    dof_map.local_variable_indices(dof_ids, mesh, scalar.first);
//
//  }

  // for conserved particles, we put a chemical potential based on the
  // prescribed density and the volume
//  for (auto&& scalar : _conservation)
//  {
//    unsigned int dof = scalar.first;
//    vector<dof_id_type> dof_ids;
//    dof_map.local_variable_indices(dof_ids, mesh, scalar.first);
//
//  }

  // for conserved particles, we put a chemical potential based on the
  // prescribed density and the volume
//  for (auto&& scalar : _conservation)
//  {
//    unsigned int dof = scalar.first;
//    vector<dof_id_type> dof_ids;
//    dof_map.local_variable_indices(dof_ids, mesh, scalar.first);
//
//  }

  // guess is checked: ok
  solution_u.close();
  poisson.update();
}


void
DriftDiffusion::do_set_to_remembered_solution(ID id)
{
  // call the default implementation
  SimulationInterface::do_set_to_remembered_solution(id);

  get_environment().prepare_for_solve();
  //build_local_scaling();

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
  }
  else if (coupling == "holes")
  {
    myopts.coupling = HCURRENT | POISSON;
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

  // setup the electrochemical potentials
  for (auto&& name : _carriers)
  {
    bool constant_qfermi = false;

    auto physit(get_options().submodels_begin("Physics"));
    ModelOptions& physopts = physit->second;

    auto itc (physopts.submodels_begin("carrier"));
    auto end_itc (physopts.submodels_end("carrier"));
    for ( ; itc != end_itc; ++itc)
    {
      if (name == itc->second.get_name())
      {
        if (itc->second.get_option("constant_qFermi", false))
          constant_qfermi = true;
      }
    }

    if (constant_qfermi)
    {
      Messages::newline();
      Messages::info("Carrier " + name +
          " will be defined with constant quasi Fermi level.");
      system.add_variable(name, libMeshEnums::FIRST,
          libMeshEnums::SCALAR, &(_carrier_region_ids[name]));
    }
    else
      system.add_variable(name, libMeshEnums::FIRST,
          &(_carrier_region_ids[name]));
  }

  // this is the last variable, so that we can use the index of the carriers also
  // as their variable index
  // the potential is assumed to be on the whole domain
  system.add_variable("potential", libMeshEnums::FIRST, &(this->get_region_ids()));

  //
  // setup conservation
  //

  auto physit(get_options().submodels_begin("Physics"));
  ModelOptions& physopts = physit->second;

  set<string> conserved_carriers;

  auto itc (physopts.submodels_begin("number_conservation"));
  auto end_itc (physopts.submodels_end("number_conservation"));

  for ( ; itc != end_itc; ++itc)
  {
    const ModelOptions& opts = itc->second;
    string name = opts.get_option("carrier", "");

    if (conserved_carriers.count(name))
      Messages::warning("Redefining number conservation for carrier " + name);

    if (name.empty())
      throw InitFailedException("You must specify a carrier for number conservation.");

    conserved_carriers.insert(name);

    unsigned int var = system.add_variable(name + "_c",
        libMeshEnums::FIRST, libMeshEnums::SCALAR);

    _conservation[var].name = name + "_c";
    _conservation[var].id = var;
    _conservation[var].options = itc->second;
    _conservation[var].carrier_vars.push_back(system.variable_number(name));
    _conservation[var].stoichiometry.push_back(1.0);

    _conservation[var].conserved_number = opts.get_option("conserved_number", 0.0);
    _conservation[var].initial_density = opts.get_option("initial_density", 0.0);

    if ((_conservation[var].conserved_number == 0.0) &&
        (_conservation[var].initial_density == 0.0))
      throw InitFailedException(string("You must specify one of 'conserved_number' or ") +
          "'initial_density' for number conservation of carrier '" + name + "'");
  }


  system.add_vector("old_sol", true, GHOSTED);
  system.add_vector("weight", true, GHOSTED);
  system.add_vector("scaling", true, GHOSTED);
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

  ostringstream name;
  name << "__DD_rstf" << dd->get_id();

  sys = &(es.add_system<RSTFSys>(name.str()));
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
  auto bdit(_boundaries.begin());
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
DriftDiffusion::RSTFSys::get_testfunction(const string& bd)
{
  libMesh::NumericVector<double>* vec = NULL;

  auto bdit(_boundaries.find(bd));
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
  auto bdit(_boundaries.begin());
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
        auto bdit(_boundaries.find(bd->get_name()));
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
    auto bdit(_boundaries.begin());
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
  data_output.write_nodal_data(this->name(), sol, solname);

}


void
DriftDiffusion::RSTFSys::build_nodal_results(vector<double>& results,
    vector<string>& legend)
{
  legend.resize(_boundaries.size());

  vector<libMesh::NumericVector<double>*> rhsides(_boundaries.size());
  auto bdit(_boundaries.begin());
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
      auto bdit(_boundaries.begin());
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


  {
    // find the region IDs for the different variables
    // this is not elegant, but it cannot be done in the constructor
    auto physit(get_options().submodels_begin("Physics"));
    ModelOptions& physopts = physit->second;

    auto itc (physopts.submodels_begin("carrier"));
    auto end_itc (physopts.submodels_end("carrier"));

    for ( ; itc != end_itc; ++itc)
    {
      string name = (itc->second).get_option("name", "");
      string regions = (itc->second).get_option("regions", get_option("regions", "all"));
      set<ID> reg_ids;
      this->get_environment().get_device().extract_physical_regions(regions, reg_ids);
      _carrier_region_ids[name].insert(reg_ids.begin(), reg_ids.end());
    }


/*
    // check excitons regions
    itc = physopts.submodels_begin("carrier");

    for ( ; itc != end_itc; ++itc)
    {
      bool ex = (itc->second).get_option("exciton", false); 
      if (ex) //the carrier is an exciton
      {
        string name = (itc->second).get_option("name", "");

        vector<string> ex_carriers;
        (itc->second).get_option("exciton_carriers", ex_carriers);
        ex_carriers.resize(2);

        vector<ID> regions;
        copy(_carrier_region_ids[name].begin(), _carrier_region_ids[name].end(), back_inserter(regions));

        for (auto exc : ex_carriers)
        {
          vector<ID> regtmp = regions;
          regions.resize(0);
          set_intersection(regtmp.begin(), regtmp.end(), 
                           _carrier_region_ids[exc].begin(), _carrier_region_ids[exc].end(), 
                           back_inserter(regions));
        }

        if (regions.size() == 0)
          throw InitFailedException("Exciton '" + name + "' and its carriers must be defined on the same regions");

      }
    }
    // end excitons
*/


  }

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
        _boundary_currents[bnd->get_name()] = 0.0;
        _voltages[bnd->get_name()] = 0.0;
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
}



void
DriftDiffusion::do_setup_solution_variables(void)
{


  // declare solution variables
  declare_solution(ElPotential, REAL, NODES, "V");
  declare_solution(ElField, VECTOR, CELL, "V/cm");

  _qFermi_base = LAST;
  _refenergy_base = _qFermi_base + _carriers.size();
  _density_base = _refenergy_base + _carriers.size();
  _mobility_base = _density_base + _carriers.size();
  _conductivity_base = _mobility_base + _carriers.size();
  _flux_base = _conductivity_base + _carriers.size();
  _curr_base = _flux_base + _carriers.size();
  _joule_base = _curr_base + _carriers.size();
  _thelpower_base = _joule_base + _carriers.size() + 1;
  _peltier_base = _thelpower_base + _carriers.size();
  _recheat_base = _peltier_base + _carriers.size() + 1; // we put them all together
  _powerflux_base = _recheat_base + 1;
  _net_rec_base = _powerflux_base + _carriers.size() + 1;
  _rec_base = _net_rec_base + _carriers.size();
  //not used now, just a reminder of the correct way to proceed adding new variables if needed
  //_next_var_base = _rec_base + _rec_models.size()*_carriers.size();

  for (unsigned int i = 0; i < _carriers.size(); ++i)
  {
    declare_solution_ext(_carriers[i] + "QFermi", _qFermi_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "eV");
    if (plot_solution("QFermi"))
      add_plot_variable(_qFermi_base + i);

    declare_solution_ext(_carriers[i] + "Density", _density_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "cm^-3");
    if (plot_solution("Density"))
          add_plot_variable(_density_base + i);

    declare_solution_ext(_carriers[i] + "RefEnergy", _refenergy_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "eV");
    if (plot_solution("RefEnergy"))
          add_plot_variable(_refenergy_base + i);

    declare_solution_ext(_carriers[i] + "Mobility", _mobility_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "cm^2/(V*s)");
    if (plot_solution("Mobility"))
          add_plot_variable(_mobility_base + i);

    declare_solution_ext(_carriers[i] + "Conductivity", _conductivity_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "S/cm");
    if (plot_solution("Conductivity"))
          add_plot_variable(_conductivity_base + i);

    declare_solution_ext(_carriers[i] + "Flux", _flux_base + i,
        SolutionDescriptor::VECTOR, SolutionDescriptor::CELL, "1/(s*cm^2)");
    if (plot_solution("Flux"))
          add_plot_variable(_flux_base + i);

    declare_solution_ext(_carriers[i] + "CurrentDensity", _curr_base + i,
        SolutionDescriptor::VECTOR, SolutionDescriptor::CELL, "A/cm^2");
    if (plot_solution("CurrentDensity"))
          add_plot_variable(_curr_base + i);

    declare_solution_ext(_carriers[i] + "ThElPower", _thelpower_base + i, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "eV/K");
    if (plot_solution("ThermoelectricPower"))
          add_plot_variable(_thelpower_base + i);

    declare_solution_ext(_carriers[i] + "Joule", _joule_base + i, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "W/cm^3");
    if (plot_solution("JouleHeat"))
          add_plot_variable(_joule_base + i);

    //declare_solution_ext(_carriers[i] + "PowerFlux", _powerflux_base + i, SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "W/cm^2");

    //declare_solution_ext(_carriers[i] + "Peltier", _peltier_base + i, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "W/cm^3");



    declare_solution_ext(_carriers[i] + "NetRecombination", _net_rec_base + i,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "1/(s*cm^3)");
    if (plot_solution("NetRecombination"))
          add_plot_variable(_net_rec_base + i);

  }

  declare_solution_ext("TotalRecombinationHeat", _recheat_base, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "W/cm^3");
  declare_solution_ext("TotalJouleHeat", _joule_base + _carriers.size(), SolutionDescriptor::REAL, SolutionDescriptor::NODES, "W/cm^3");
  if (plot_solution("JouleHeat"))
    add_plot_variable(_joule_base + _carriers.size());
  //declare_solution_ext("TotalPowerFlux", _powerflux_base + _carriers.size(), SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "W/cm^2");
  //declare_solution_ext("TotalPeltierHeat", _peltier_base + _carriers.size(), SolutionDescriptor::REAL, SolutionDescriptor::NODES, "W/cm^3");

  unsigned int n_rec = 0;
  for (auto& rec : _rec_models)
  {
    for (auto id : rec.second)
    {
      declare_solution_ext(_carriers[id] + "Recombination_" + rec.first, _rec_base + n_rec*_carriers.size() + id,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES, "1/(s*cm^3)");

      if (plot_solution("NetRecombination"))
        add_plot_variable(_rec_base + n_rec*_carriers.size() + id);
    }
    n_rec++;
  }

  /*
  // the correct number of components will be inserted afterwards
  declare_solution(ElectronBands, NTUPLE, CELL, "eV", 1);
  declare_solution(HoleBands, NTUPLE, CELL, "eV", 1);
  if (plot_solution("BandEdges"))
  {
    add_plot_variable(ElectronBands);
    add_plot_variable(HoleBands);
  }
  */

  declare_solution(Polarization, VECTOR, CELL, "C/m^2");

  declare_solution(TotalCurrentDensity, VECTOR, CELL, "A/cm^2");
  if (plot_solution("CurrentDensity"))
    add_plot_variable(TotalCurrentDensity);

  declare_solution(IonizedDonors, REAL, NODES, "cm^-3");
  declare_solution(IonizedAcceptors, REAL, NODES, "cm^-3");

  /*
  declare_solution(IonizedElectronTraps, REAL, NODES, "cm^-3");
  declare_solution(IonizedHoleTraps, REAL, NODES, "cm^-3");



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
  */

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

  const unsigned int n_vars = _carriers.size();
  // the carriers
  vector<unsigned int> q_var_num(n_vars);
  for (unsigned int i = 0; i < n_vars; ++i)
    q_var_num[i] = system->variable_number(_carriers[i]);


  libMesh::FEType fe_type = system->variable_type(u_var);
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<vector<unsigned int>> dof_indices_q(n_vars);

  map<unsigned int, vector<unsigned int>> dof_indices_cons;


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

  const unsigned int n_dofs = dof_indices_u.size();

  // to handle constant qFermi case transparently
  map<unsigned int, vector<unsigned int> > id_map;

  // these are the variables present in the element
  set<unsigned int> q_var;
  for (unsigned int i = 0; i < n_vars; ++i)
  {
    dof_map.dof_indices(elem, dof_indices_q[i], q_var_num[i]);
    if (dof_indices_q[i].size() > 0)
      q_var.insert(i);

    id_map[i].resize(n_dofs);
    if (dof_indices_q[i].size() == 1)
      fill(id_map[i].begin(), id_map[i].end(), 0);
    else
      iota(id_map[i].begin(), id_map[i].end(), 0);
  }

  for (auto&& scalar : _conservation)
  {
    unsigned int dof = scalar.first;
    dof_map.dof_indices(elem, dof_indices_cons[dof], dof);
  }

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
  libMesh::RealGradient curr_tot(0);
  vector<RealGradient> flux(_carriers.size(), 0);
  vector<RealGradient> curr(_carriers.size(), 0);

  libMesh::RealGradient jn(0);
  libMesh::RealGradient jp(0);
  libMesh::RealGradient el_field(0);
  libMesh::RealVectorValue polariz(0);
  double el_pot;

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    vector<double> qf(n_vars, 0.0);
    vector<double> pot_erg(n_vars, 0.0);
    double oldu  = 0.0;
    vector<double> oldqf(n_vars, 0.0);
    double T  = 0.0;
    vector<double> thel_pow(n_vars, 0.0);

    libMesh::RealGradient e_field(0);
    vector<libMesh::RealGradient> grad_qf_loc(n_vars, 0);
    libMesh::RealGradient grad_T_loc(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices_u[i]);
      oldu += phi[i][n] * oldsolution(dof_indices_u[i]);
      e_field += dphi[i][n] * solution(dof_indices_u[i]);

      for (auto& v : q_var)
      {
        unsigned int ii = id_map[v][i];
        qf[v] += phi[i][n] * solution(dof_indices_q[v][ii]);
        oldqf[v] += phi[i][n] * oldsolution(dof_indices_q[v][ii]);
        grad_qf_loc[v] += dphi[i][n] * solution(dof_indices_q[v][ii]);
      }


      grad_T_loc += dphi[i][n] * T_nodes[i];

      T +=  phi[i][n] * T_nodes[i];

    }

    // scale the potential back
    u *= phi0;
    oldu *= phi0;
    e_field *= -phi0;

    for (auto&& cons : _conservation)
    {
      for (auto&& var : cons.second.carrier_vars)
      {
        // only if the carrier is present in the element
        if (q_var.count(var))
        {
          qf[var] += solution(dof_indices_cons[cons.first][0]);
          oldqf[var] += oldsolution(dof_indices_cons[cons.first][0]);
        }
      }
    }

    for (auto& v : q_var)
    {
      qf[v] *= phi0;
      oldqf[v] *= phi0;
      grad_qf_loc[v] *= phi0;
      pot_erg[v] = -u * std::fabs(sc->get_carrier_properties(v)->get_charge());
    }

    el_pot = u;


    sc->set_coordinates(real_pts[n]);

    sc->set_el_potential(u);
    sc->set_old_el_potential(oldu);
    sc->set_electric_field(e_field);

    for (auto& v : q_var)
    {
      sc->set_fermi_potential(v, qf[v]);
      sc->set_old_fermi_potential(v, oldqf[v]);
      sc->set_grad_fermi(v, grad_qf_loc[v]);
    }

    sc->calculate_densities();
    sc->calculate_mobilities();

    // NOTE: from the definition of the (electro)chemical potential, uncharged
    // carriers are consistent with negatively charged carriers in terms of flux
    for (auto& v : q_var)
    {
      double q = sc->get_carrier_properties(v)->get_charge();
      double sign = sc->get_carrier_properties(v)->get_charge_sign();
      double sigma = sc->get_q_conductivity(v);
      thel_pow[v] = sc->get_carrier_properties(v)->get_thermoelectric_power();
      RealGradient flux_loc = -sigma * (sign * grad_qf_loc[v] + thel_pow[v] * grad_T_loc);
      flux[v] += flux_loc;
      curr[v] += Constants::e * q  * flux_loc;
      curr_tot += Constants::e * q * flux_loc;
    }


    el_field += e_field;
    polariz += sc->get_total_polarization();


    if (values.count(ElPotential))
      values[ElPotential][n] = u - _reference_potential;

    for (unsigned int v = 0; v < n_vars; ++v)
    {
      if (values.count(_qFermi_base + v))
        values[_qFermi_base + v][n] = -qf[v];

      if (values.count(_density_base + v))
        values[_density_base + v][n] = sc->get_q_density(v);

      if (values.count(_refenergy_base + v))
        values[_refenergy_base + v][n] = sc->get_carrier_band_edge(v) + pot_erg[v];

      if (values.count(_mobility_base + v))
        values[_mobility_base + v][n] = sc->get_q_mobility(v);

      if (values.count(_conductivity_base + v))
      {
        double q = 0;
        double sign = 1;
        if (q_var.count(v))
        {
          q = sc->get_carrier_properties(v)->get_charge();
          sign = sc->get_carrier_properties(v)->get_charge_sign();
        }
        values[_conductivity_base + v][n] = Constants::e * q * sign * sc->get_q_conductivity(v);
      }

      if (values.count(_thelpower_base + v))
        values[_thelpower_base + v][n] = thel_pow[v];

    }

    /*
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

    */



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


    if (values.count(_joule_base + n_vars))
      values[_joule_base + n_vars][n] = 0.0;

    for (unsigned int v = 0; v < n_vars; ++v)
    {
      if (sc->get_carrier_properties(v) != nullptr)
      {
        double q = sc->get_carrier_properties(v)->get_charge();
        double sign = sc->get_carrier_properties(v)->get_charge_sign();
        double joule_loc = -sign * Constants::e * (flux[v] * grad_qf_loc[v]);
        if (values.count(_joule_base + v))
        {
          values[_joule_base + v][n] = joule_loc;
        }

        if (values.count(_joule_base + n_vars))
        {
          values[_joule_base + n_vars][n] += joule_loc;
        }
      }
      else
      {
        if (values.count(_joule_base + v))
          values[_joule_base + v][n] = 0;
      }
    }

    /*
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


    */

    bool need_recomb = false;

    if (values.count(_recheat_base))
    {
      values[_recheat_base][n] = 0;
      need_recomb = true;
    }

    for (unsigned int v = 0; v < n_vars; ++v)
    {
      if (values.count(_net_rec_base + v))
        need_recomb = true;
    }

    if (need_recomb)
      sc->calculate_net_recombination_rates();

    for (unsigned int v = 0; v < n_vars; ++v)
    {
      if (values.count(_net_rec_base + v))
        values[_net_rec_base + v][n] = sc->get_net_q_recombination_rate(v);
    }

    unsigned int n_rec = 0;
    for (auto& rec : _rec_models)
    {
      for (auto id : rec.second)
      {
        if (values.count(_rec_base + n_rec*_carriers.size() + id))
        {
          values[_rec_base + n_rec*_carriers.size() + id][n] = 0.0;

          DriftDiffusionProperties::RecombinationModelIterator recit(sc->recombination_models_begin());
          DriftDiffusionProperties::RecombinationModelIterator recend(sc->recombination_models_end());

          for ( ; recit != recend; ++recit)
          {
            if (rec.first == (recit->second)->get_plot_name())
            {
              vector<double> R(_carriers.size(), 0.0);
              vector<vector<double>> dR( _carriers.size(), vector<double>(_carriers.size() + 1, 0.0) );

              (recit->second)->get_net_rate_and_derivatives(R, dR);
              values[_rec_base + n_rec*_carriers.size() + id][n] += R[id];
            }
          }
        }
      }
      n_rec++;
    }


    if (values.count(_recheat_base))
    {
      DriftDiffusionProperties::RecombinationModelIterator recit(sc->recombination_models_begin());
      DriftDiffusionProperties::RecombinationModelIterator recend(sc->recombination_models_end());

      for ( ; recit != recend; ++recit)
      {
        //if (! (recit->second)->is_radiative())
        {
          vector<double> R(_carriers.size(), 0.0);
          vector<vector<double>> dR( _carriers.size(), vector<double>(_carriers.size() + 1, 0.0) );

          (recit->second)->get_net_rate_and_derivatives(R, dR);

          for (auto& v : q_var)
          {
            double sign = sc->get_carrier_properties(v)->get_charge_sign();
            double P = thel_pow[v];
            double H = Constants::e * R[v] * (sign * qf[v] + T * P);
            values[_recheat_base][n] += H;


            if ((recit->second)->is_radiative())
            {
              double energy = 0;
              values[_recheat_base][n] += sign * Constants::e * R[v] * sc->get_carrier_band_edge(v);
            }
          }

        }
      }
    }


    /*
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
    */

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

  if (values.count(TotalCurrentDensity))
  {
    values[TotalCurrentDensity][0] = curr_tot(0) / np;
    values[TotalCurrentDensity][1] = curr_tot(1) / np;
    values[TotalCurrentDensity][2] = curr_tot(2) / np;
  }

  for (unsigned int v = 0; v < n_vars; ++v)
  {
    if (values.count(_flux_base + v))
    {
      values[_flux_base + v][0] = flux[v](0) / np;
      values[_flux_base + v][1] = flux[v](1) / np;
      values[_flux_base + v][2] = flux[v](2) / np;
    }

  }

  for (unsigned int v = 0; v < n_vars; ++v)
  {
    if (values.count(_curr_base + v))
    {
      values[_curr_base + v][0] = curr[v](0) / np;
      values[_curr_base + v][1] = curr[v](1) / np;
      values[_curr_base + v][2] = curr[v](2) / np;
    }

  }
/*
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
*/
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
        Xn.reposition(n_dofs, n_dofs);
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

          _boundary_eqfermi[bd->get_name()] += J * en * phi0;
          _boundary_hqfermi[bd->get_name()] += J * ep * phi0;
          boundary_area[bd->get_name()] += J;
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
      reffermi_e[sit->first] += _boundary_eqfermi[(*bdit)->get_name()] / count;
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
  map<string, libMesh::NumericVector<double>*> rstf;

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

  const unsigned int n_vars = _carriers.size();
  // the carriers
  vector<unsigned int> qf_var_num(n_vars);
  for (unsigned int i = 0; i < n_vars; ++i)
    qf_var_num[i] = system->variable_number(_carriers[i]);


  // all have the same type
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


  // dof indices of electrostatic potential
  vector<unsigned int> dof_indices_u;
  // dof indices of all QF potentials
  vector<vector<unsigned int>> dof_indices_qf(n_vars);

  // dof indices for the RSTF functions
  vector<unsigned int> dof_indices_rstf;


  MeshBase::const_element_iterator el(this->active_local_elements_begin());
  const MeshBase::const_element_iterator end_el(this->active_local_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    //const Elem* top_parent = elem->top_parent();


    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);

    // these are the variables actually present in the element
    set<unsigned int> qf_vars;

    // a remapping for local dof ids, to allow for constant qFermi variables
    // this will allow to use the normal assembly functions, but with constant
    // potential
    map<unsigned int, vector<unsigned int> > id_map;

    for (unsigned int i = 0; i < n_vars; ++i)
    {
      dof_map.dof_indices(elem, dof_indices_qf[i], qf_var_num[i]);
      if (dof_indices_qf[i].size() > 0)
        qf_vars.insert(i);

      id_map[i].resize(dof_indices_u.size());

      if (dof_indices_qf[i].size() == 1)
        fill(id_map[i].begin(), id_map[i].end(), 0);
      else
        iota(id_map[i].begin(), id_map[i].end(), 0);

    }

    dof_map_rstf.dof_indices(elem, dof_indices_rstf, 0);


    DDBulkModel* sc =
        get_bulk_model<DDBulkModel>(elem);


    assert(sc != NULL);


    fe->reinit(elem);

    sc->reinit(elem);

    //Get the temperature given the element
    vector<double> T_nodes =  sc->get_temperature_at_nodes();


    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      // this is the same for all
      unsigned int n_dofs = dof_indices_u.size();

      Real u  = 0.0;
      libMesh::RealGradient e_field(0);

      vector<Real> qf(n_vars, 0.0);
      vector<RealGradient> grad_qf(n_vars, 0);

      RealGradient dT(0);

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        e_field -= dphi[i][qp] * solution(dof_indices_u[i]);

        for (auto& var : qf_vars)
        {
          unsigned int ii = id_map[var][i];
          qf[var] += phi0 * phi[i][qp] * solution(dof_indices_qf[var][ii]);

          grad_qf[var] += dphi[i][qp] * phi0 * solution(dof_indices_qf[var][ii]);
        }

        dT += dphi[i][qp] * T_nodes[i];
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);


      sc->set_el_potential(phi0 * u);
      sc->set_electric_field(phi0 * e_field);

      for (auto& var : qf_vars)
      {
        sc->set_fermi_potential(var, qf[var]);
        sc->set_grad_fermi(var, grad_qf[var]);
      }


      sc->calculate_densities();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      double net_rate = 0.0;
      vector<double> sigma(n_vars, 0.0);
      for (auto& var : qf_vars)
      {
        double chrg = sc->get_carrier_properties(var)->get_charge();
        sigma[var] = Constants::e * chrg * sc->get_q_conductivity(var);

        double Rn = sc->get_net_q_recombination_rate(var);
        net_rate -= chrg * Rn;
      }

      net_rate *= JxW[qp] * Constants::e;

      // loop over all recombination models and add the ones that model
      // tunneling for the current contact
      // this is needed, because the current calculated afterwards does not
      // include direct tunneling.
      /*
       * TODO: this has to be adjusted

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
      */

      vector<RealGradient> curr(n_vars);
      for (auto& var : qf_vars)
      {
        double sign = sc->get_carrier_properties(var)->get_charge_sign();
        double P = sc->get_carrier_properties(var)->get_thermoelectric_power();
        curr[var] = JxW[qp] * sigma[var] * (sign * grad_qf[var] + P * dT);
      }

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        // do this for each contact
        auto rstf_it(rstf.begin());
        for ( ; rstf_it != rstf.end(); ++rstf_it)
        {
          const libMesh::NumericVector<double>& sol = *rstf_it->second;

          for (auto& var : qf_vars)
          {
            double value = (curr[var] * dphi[n][qp] +
                net_rate * phi[n][qp]) * sol(dof_indices_rstf[n]);
            _boundary_currents[rstf_it->first] -= value;
          }

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

  const unsigned int n_vars = _carriers.size();
  // the carriers
  vector<unsigned int> qf_var_num(n_vars);
  for (unsigned int i = 0; i < n_vars; ++i)
    qf_var_num[i] = system->variable_number(_carriers[i]);

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
  // dof indices of all QF potentials
  vector<vector<unsigned int>> dof_indices_qf(n_vars);


  BoundaryElementMap::iterator el(env.boundary_elements_begin());
  BoundaryElementMap::iterator end_el(env.boundary_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    //const Elem* top_parent = elem->top_parent();

    const Boundary* boundary = el.get_boundary();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_u, u_var);

    // these are the variables actually present in the element
    set<unsigned int> qf_vars;

    // a remapping for local dof ids, to allow for constant qFermi variables
    // this will allow to use the normal assembly functions, but with constant
    // potential
    map<unsigned int, vector<unsigned int> > id_map;

    for (unsigned int i = 0; i < n_vars; ++i)
    {
      dof_map.dof_indices(elem, dof_indices_qf[i], qf_var_num[i]);
      if (dof_indices_qf[i].size() > 0)
        qf_vars.insert(i);

      id_map[i].resize(dof_indices_u.size());

      if (dof_indices_qf[i].size() == 1)
        fill(id_map[i].begin(), id_map[i].end(), 0);
      else
        iota(id_map[i].begin(), id_map[i].end(), 0);

    }

    DDBulkModel* sc =
        get_bulk_model<DDBulkModel>(elem);

    assert(sc != NULL);


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
      libMesh::RealGradient e_field(0);

      vector<Real> qf(n_vars, 0.0);
      vector<RealGradient> grad_qf(n_vars, 0);

      libMesh::RealGradient dT(0);

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        e_field -= dphi[i][qp] * solution(dof_indices_u[i]);


        for (auto& var : qf_vars)
        {
          unsigned int ii = id_map[var][i];
          qf[var] += phi0 * phi[i][qp] * solution(dof_indices_qf[var][ii]);

          grad_qf[var] += dphi[i][qp] * phi0 * solution(dof_indices_qf[var][ii]);
        }

        dT += dphi[i][qp] * T_nodes[i];
      }

      // prepare for calculating local properties
      //sc->set_coordinates(elem->centroid());  ????? 2012-08-31
      sc->set_coordinates(q_point[qp]);

      sc->set_el_potential(phi0 * u);
      sc->set_electric_field(phi0 * e_field);

      for (auto& var : qf_vars)
      {
        sc->set_fermi_potential(var, qf[var]);
        sc->set_grad_fermi(var, grad_qf[var]);
      }

      sc->calculate_densities();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      //Get the thermoelectric power
      /*
      sc->compute_thermoelectric_powers();
      double Pn =  sc->get_electron_thermoelectric_power() / phi0;
      double Pp =  sc->get_hole_thermoelectric_power() / phi0;
      */

      double net_rate = 0.0;
      vector<RealGradient> curr(n_vars, 0.0);
      for (auto& var : qf_vars)
      {
        double chrg = sc->get_carrier_properties(var)->get_charge();
        double sign = sc->get_carrier_properties(var)->get_charge_sign();

        double sigma = Constants::e * chrg * sc->get_q_conductivity(var);
        curr[var] = JxW[qp] * sigma * sign * grad_qf[var];

        double Rn = sc->get_net_q_recombination_rate(var);
        net_rate -= chrg * Rn;
      }

      net_rate *= JxW[qp] * Constants::e;

      //libMesh::RealGradient je(JxW[qp] * phi0 * (sigma_e * (dEfn + Pn * dT)));
      //libMesh::RealGradient jh(JxW[qp] * phi0 * (sigma_h * (dEfp + Pp * dT)));

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        for (auto& var : qf_vars)
        {
          _boundary_currents[boundary->get_name()] -= (curr[var] * dphi[n][qp] +
              net_rate * phi[n][qp]) * weight[n];
        }
      }

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
        fe_currents[bnd->get_name()] = 0.0;
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

        fe_currents[boundary->get_name()] += current;
      }
    } // end loop over elem sides
  } // end loop over elements


  Messages m;
  m.info("Field emission currents:");
  m.indent();
  for (ContactData::iterator it(fe_currents.begin()); it != fe_currents.end(); ++it)
  {
    ostringstream os;
    os << (*it).first << ": " << (*it).second << "" << endl;
    m.info(os.str());
  }
}



void
DriftDiffusion::calculate_currents_surfint(void)
{
  Messages::warning("Current by surface integration is not implemented.\n");

/*
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

*/
}


void
DriftDiffusion::calculate_surface_recombination(void)
{
  //Messages::warning("Current by surface integration is not implemented.\n");

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  // need the ghosted solution vector
  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const MeshBase& mesh = system->get_mesh();
  const Device& device = *(_device);
  const SimulationEnvironment& env = get_environment();

  const libMesh::DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  const double phi0 = get_scaling().get_potential_scaling();


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");

  const unsigned int n_vars = _carriers.size();
  // the carriers
  vector<unsigned int> qf_var_num(n_vars);
  for (unsigned int i = 0; i < n_vars; ++i)
    qf_var_num[i] = system->variable_number(_carriers[i]);

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
  // dof indices of all QF potentials
  vector<vector<unsigned int>> dof_indices_qf(n_vars);

  // the total recombination current
  double current = 0.0;

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

    // these are the variables actually present in the element
    set<unsigned int> qf_vars;
    map<unsigned int, vector<unsigned int> > id_map;

    for (unsigned int i = 0; i < n_vars; ++i)
    {
      dof_map.dof_indices(elem, dof_indices_qf[i], qf_var_num[i]);
      if (dof_indices_qf[i].size() > 0)
        qf_vars.insert(i);

      id_map[i].resize(dof_indices_u.size());

      if (dof_indices_qf[i].size() == 1)
        fill(id_map[i].begin(), id_map[i].end(), 0);
      else
        iota(id_map[i].begin(), id_map[i].end(), 0);

    }

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
          libMesh::RealGradient e_field(0);

          vector<Real> qf(n_vars, 0.0);
          vector<RealGradient> grad_qf(n_vars, 0);

          for (unsigned int i = 0; i < phi_size; i++)
          {
            u  += phi[i][qp] * solution(dof_indices_u[i]);
            e_field += dphi[i][qp] * solution(dof_indices_u[i]);

            double tmp = dphi[i][qp] * face_normals[qp];
            for (auto& var : qf_vars)
            {
              unsigned int ii = id_map[var][i];
              qf[var] += phi[i][qp] * phi0 * solution(dof_indices_qf[var][ii]);

              grad_qf[var] += tmp * dphi[i][qp] * phi0 * solution(dof_indices_qf[var][ii]);
            }

          }

          // prepare for calculating local properties
          sm->set_coordinates(q_point[qp]);
          sm->set_el_potential(phi0 * u);
          sm->set_electric_field(-phi0 * e_field);

          for (auto& var : qf_vars)
          {
            sm->set_fermi_potential(var, qf[var]);
            sm->set_grad_fermi(var, grad_qf[var]);
          }

          sm->set_face_normal(face_normals[qp]);
          sm->compute();

          //vector<double> value_v(_carriers.size() + 1, 0.0);

          const vector<double>& coeff_a = sm->get_a();
          const vector<double>& coeff_g = sm->get_g();

          //double tmp = (coeff_g[0] - coeff_a[0] * u * phi0);
          //if (sm->get_type(0) != DDInterfaceModel::DIRICHLET)
          //  value_u += J * tmp / (x0 * C0);
          //else
          //  value_u = coeff_g[0] / phi0;

          for (auto&& var : qf_vars)
          {
            if (var != u_var)
            {
              double value = 0.0;
              if (sm->get_type(var) != DDInterfaceModel::DIRICHLET)
                value = (coeff_g[var] - coeff_a[var] * qf[var] * phi0);

              current += JxW[qp] * 0.5 * value;
            }
          }

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

  const libMesh::NumericVector<Number>& solution = system.get_solution_vector();
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
  double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;
  // scaling for recombination rates
  double R0 = C0 / scaling.get_time_scaling();

  double C0_q = C0;
  C0_q = 1.0;

  //cout<<"l2 = "<<l2<<" C0 = "<<C0<<" x0 = "<<x0<<" mu0 = "<<mu0<<" C0_q = "<<C0_q<<endl;

  const unsigned int u_var = system.variable_number("potential");
  // fermi potential variable numbers are defined within element loop
  // since different variables can ben set in different regions

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
  // since different variables can ben set in different regions

  auto it = this->active_local_elements_begin();
  const auto end = this->active_local_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    DDBulkModel* sc = get_bulk_model<DDBulkModel>(elem);
    assert(sc != NULL);

    sc->reinit(elem);
    fe->reinit(elem);

    //Get variables for fermi potentials
    set<unsigned int> q_var;
    for (auto&& cp : sc->get_carrier_properties())
      q_var.insert(cp.first);

    //Get dof indices
    map<unsigned int, vector<unsigned int>> dof_indices_q;

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    dof_indices.resize(0);
    for (auto var : q_var)
    {
      dof_map.dof_indices(elem, dof_indices_q[var], var);
      dof_indices.insert(dof_indices.end(),
          dof_indices_q[var].begin(), dof_indices_q[var].end());
    }
    dof_indices.insert(dof_indices.end(),
        dof_indices_u.begin(), dof_indices_u.end());

    local_scaling.resize(dof_indices.size());

    unsigned int n_dofs = dof_indices_u.size();

    // a remapping for local dof ids, to allow for constant qFermi variables
    // this will allow to use the normal assembly functions, but with constant
    // potential
    map<unsigned int, vector<unsigned int> > id_map;
    for (auto var : q_var)
    {
      id_map[var].resize(n_dofs);
      if (dof_indices_q[var].size() == 1)
        fill(id_map[var].begin(), id_map[var].end(), 0);
      else
        iota(id_map[var].begin(), id_map[var].end(), 0);
    }
    id_map[u_var].resize(n_dofs);
    iota(id_map[u_var].begin(), id_map[u_var].end(), 0);

    //Define scale vector for fermi potentials
    map<unsigned int, libMesh::DenseSubVector<Number>> scaleq;
    unsigned int n_var = 0;
    for (auto var : q_var)
    {
      unsigned int var_dofs = dof_indices_q[var].size();
      libMesh::DenseSubVector<Number> scale_tmp(local_scaling);
      scale_tmp.reposition(n_var, var_dofs);
      scaleq.insert( make_pair(var, scale_tmp) );
      n_var += var_dofs;
    }
    scaleu.reposition(n_var, n_dofs);


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
        q.insert( make_pair(var, 0.0) );

      libMesh::RealGradient e_field(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        for (auto var : q_var)
          q[var] += phi[i][qp] * solution(dof_indices_q[var][id_map[var][i]]);

        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }

      //cout<<"u = "<<u;
      //for (auto var : q_var)
        //cout<<" q"<<var<<" = "<<q[var];

      //cout<<endl;

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_el_potential(phi0 * u);
      for (auto var : q_var)
      {
        sc->set_fermi_potential(var, phi0 * q[var]);
        sc->set_grad_fermi(var, libMesh::RealGradient(0.0));
      }
      sc->set_electric_field(phi0 / x0 * e_field);

      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      const libMesh::RealTensor& permittivity = sc->get_relative_permittivity();

      double l2_eps = JxW[qp] * l2;

      map<unsigned int, double> sigma, ddens_dphi;
      map<unsigned int, map<unsigned int, long double>> dR;

      for (auto var : q_var)
      {
        sigma.insert( make_pair(var, JxW[qp] * sc->get_q_conductivity(var) / mu0 / C0_q) );
        ddens_dphi.insert( make_pair(var, sc->get_q_density_derivative(var)) );
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

      double drho = 0.0;

      for (auto&& var : q_var)
      {
        if (var != u_var)
        {
          double drho_v = sc->get_charge_density_derivative(var) * phi0 / C0;
          drho += drho_v;

          for (auto&& varj : q_var)
          {
            long double dR_dEf_tmp = sc->get_net_q_recombination_rate_derivatives(var)[varj];
            long double dR_tmp = dR_dEf_tmp * phi0 / R0;

            dR[var][varj] = dR_tmp;
          }

        }
      }
      drho *= -JxW[qp];

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        double phi_i_x_phi_j = JxW[qp] * phi[i][qp] * phi[i][qp];

        for (auto&& scale : scaleq)
        {
          const char ct = sc->get_carrier_properties(scale.first)->get_carrier_type();
          double sign = (ct == 'e') ? 1.0 : -1.0;
          scale.second(id_map[scale.first][i]) += sigma[scale.first] * (dphi[i][qp] * dphi[i][qp])
                             -sign * dR[scale.first][scale.first] * phi_i_x_phi_j;
        }

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
            for (auto&& dofs : dof_indices_q)
            {
              unsigned int ii = id_map[dofs.first][i];
              if (sm->get_type(dofs.first) == DDInterfaceModel::DIRICHLET)
              dirichlet_dofs.insert(dofs.second[ii]);
            }

            if (sm->get_type(u_var) == DDInterfaceModel::DIRICHLET)
              dirichlet_dofs.insert(dof_indices_u[i]);

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
    //for (size_t i = 0; i< local_scaling.size(); i++ )
      //cout<<"ls["<<i<<"] = "<<local_scaling(i)<<endl;

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
  loc_scaling.localize(loc_scaling, system.get_dof_map().get_send_list());
  //loc_scaling.print_matlab("scaling.m");


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
        if (tokens[0] == it->first)
        {
          const Boundary* bnd = this->get_environment().get_boundary(it->first);
          double curr = it->second * bnd->get_area_factor();
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
    case (CURRENTS):
      _this->do_assembly<CURRENTS>(x, residual, jacobian);
      break;
    case (POISSON):
      _this->do_assembly<POISSON>(x, residual, jacobian);
      break;
    default:
      _this->do_assembly<FULLYCOUPLED>(x, residual, jacobian);
      break;
  }

}



void
DriftDiffusion::do_check_nonlinear_step(
        libMesh::NumericVector<Number>& dx)
{
  if (!get_option("limit_step", false)) return;
  // references for nicer code
  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const unsigned int dim = mesh.mesh_dimension();

  //const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  const libMesh::DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");

  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  // dof indices vectors for fermi potentials are defined within element loop
  // since different variables can ben set in different regions



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

    // Get variables for fermi potentials
    set<ID> q_var;

    // we will loop only over carriers that are present in this element
    q_var.insert(u_var);
    for (auto&& cp : sc->get_carrier_properties())
      q_var.insert(cp.first);


    //Get dof indices
    map<unsigned int, vector<unsigned int>> dof_indices_var;

    dof_indices.clear();

    for (auto var : q_var)
    {
      //insert variable number and dof indices for the element
      dof_map.dof_indices(elem, dof_indices_var[var], var);
      dof_indices.insert(dof_indices.end(),
          dof_indices_var[var].begin(), dof_indices_var[var].end());
    }

    for (auto&& scalar : _conservation)
    {
      unsigned int dof = scalar.first;
      dof_map.dof_indices(elem, dof_indices_var[dof], dof);
      dof_indices.insert(dof_indices.end(),
          dof_indices_var[dof].begin(), dof_indices_var[dof].end());
    }

    // dof_indices now contains all DOFs, also scalar ones

    // they have all the same number of DOFs
    unsigned int n_dofs     = dof_indices_var[u_var].size();
    unsigned int n_dofs_tot = dof_indices.size();

    for (unsigned int i = 0; i < n_dofs; ++i)
    {
      double dmax = 0;
      for (auto& cp : sc->get_carrier_properties())
      {
        ID var = cp.first;
        double chrg = sc->get_carrier_properties(var)->get_charge();

        if (chrg != 0) chrg = fabs(chrg);
        double chempot = fabs(dx(dof_indices_var[var][i]) -
          chrg * dx(dof_indices_var[u_var][i]));

        if (chempot > dmax) dmax = chempot;
      }

      double maxstep = 10;
      if (dmax > maxstep)
      {
        for (auto var : q_var)
        {
          double value = maxstep * dx.el(dof_indices_var[var][i]) / dmax;
          dx.set(dof_indices_var[var][i], value);
        }

      }
    }
  }
  dx.close();
}



template <int coupling>
void
DriftDiffusion::do_assembly(const libMesh::NumericVector<Number>& x,
    libMesh::NumericVector<Number>* residual,
    libMesh::SparseMatrix<Number>* jacobian)
{

  START_LOG(get_name() + ": Matrix assembly", "");

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
  if (jacobian != NULL)
    build_local_scaling();
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
  //const double C0 = (_do_local_scaling) ? 1.0 : scaling.get_density_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  // x 1e4 because we calculate in cm, but P comes in C/m^2
  const double P0 = (Constants::e * x0 * C0) * 1e4;

  double C0_q = C0;
  C0_q = 1.0;

  // scaling for recombination rates
  double R0 = C0_q / scaling.get_time_scaling();

  const libMesh::DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
 
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

    // Get variables for fermi potentials
    set<ID> q_var;

    // we will loop only over carriers that are present in this element
    q_var.insert(u_var);
    for (auto&& cp : sc->get_carrier_properties())
      q_var.insert(cp.first);



    //Get dof indices
    map<unsigned int, vector<unsigned int>> dof_indices_var;

    dof_indices.clear();

    for (auto var : q_var)
    {
      //insert variable number and dof indices for the element
      dof_map.dof_indices(elem, dof_indices_var[var], var);
      dof_indices.insert(dof_indices.end(),
          dof_indices_var[var].begin(), dof_indices_var[var].end());
    }

    // they have all the same number of DOFs
    unsigned int n_dofs     = dof_indices_var[u_var].size();

    // a remapping for local dof ids, to allow for constant qFermi variables
    // this will allow to use the normal assembly functions, but with constant
    // potential
    map<unsigned int, vector<unsigned int> > id_map;
    for (auto var : q_var)
    {
      id_map[var].resize(n_dofs);
      if (dof_indices_var[var].size() == 1)
        fill(id_map[var].begin(), id_map[var].end(), 0);
      else
        iota(id_map[var].begin(), id_map[var].end(), 0);
    }


    for (auto&& scalar : _conservation)
    {
      unsigned int dof = scalar.first;
      dof_map.dof_indices(elem, dof_indices_var[dof], dof);
      dof_indices.insert(dof_indices.end(),
          dof_indices_var[dof].begin(), dof_indices_var[dof].end());
    }

    // dof_indices now contains all DOFs, also scalar ones

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
    map<unsigned int, map<unsigned int, DenseSubMatrix<Number>>> Kvv;
    map<unsigned int, DenseSubVector<Number>> Fv,
                                              Xv,
                                              oldXv,
                                              scalev;

    // Reposition the submatrices according to this scheme:
    //
    //        | K11 K12 ... K1u |        | F1 |
    //   Ke = | K21 K22 ... ... |;  Fe = | F2 |
    //        | ... ... ... ... |        | .. |
    //        | Ku1 Ku2 ... Kuu |        | Fu |
    //
    // Note: conservation variables are added at the end

    unsigned int n_var = 0;
    for (auto&& var : q_var)
    {

      Fv.insert(make_pair(var, DenseSubVector<Real>(Fe)));
      Xv.insert(make_pair(var, DenseSubVector<Real>(X)));
      oldXv.insert(make_pair(var, DenseSubVector<Real>(oldX)));
      scalev.insert(make_pair(var, DenseSubVector<Real>(local_scaling)));

      unsigned int var_dofs = dof_indices_var[var].size();

      Fv.at(var).reposition(n_var, var_dofs);
      Xv.at(var).reposition(n_var, var_dofs);
      oldXv.at(var).reposition(n_var, var_dofs);
      scalev.at(var).reposition(n_var, var_dofs);

      unsigned int n_varj = 0;
      for (auto&& varj : q_var)
      {
        unsigned int varj_dofs = dof_indices_var[varj].size();

        Kvv[var].insert(make_pair(varj, DenseSubMatrix<Real>(Ke)));
        Kvv[var].at(varj).reposition(n_var, n_varj, var_dofs, varj_dofs);
        n_varj += varj_dofs;
      }
      n_var += var_dofs;
    }

    unsigned int c_var = 0;
    for (auto&& cons : _conservation)
    {
      unsigned int var = cons.first;

      Fv.insert(make_pair(var, DenseSubVector<Real>(Fe)));
      Xv.insert(make_pair(var, DenseSubVector<Real>(X)));
      oldXv.insert(make_pair(var, DenseSubVector<Real>(oldX)));
      //scalev.insert(make_pair(var, DenseSubVector<Real>(local_scaling)));

      Fv.at(var).reposition(n_var + c_var, 1);
      Xv.at(var).reposition(n_var + c_var, 1);
      oldXv.at(var).reposition(n_var + c_var, 1);
      //scalev.at(var).reposition(n_var*n_dofs, 1);

      unsigned int n_varj = 0;
      for (auto&& varj : q_var)
      {
        unsigned int varj_dofs = dof_indices_var[varj].size();

        if ((find(cons.second.carrier_vars.begin(),
            cons.second.carrier_vars.end(), varj) !=
                cons.second.carrier_vars.end()) ||
            (varj == u_var))
        {
          Kvv[var].insert(make_pair(varj, DenseSubMatrix<Real>(Ke)));
          Kvv[var].at(varj).reposition(n_var + c_var,
                                       n_varj, 1, varj_dofs);
        }

        Kvv[varj].insert(make_pair(var, DenseSubMatrix<Real>(Ke)));
        Kvv[varj].at(var).reposition(n_varj,
                                     n_var + c_var, varj_dofs, 1);
        n_varj += varj_dofs;
      }

      Kvv[var].insert(make_pair(var, DenseSubMatrix<Real>(Ke)));
      Kvv[var].at(var).reposition(n_var + c_var,
                                  n_var + c_var, 1, 1);
      c_var++;

    }


    // Get the temperature given the element
    vector<double> T_nodes = sc->get_temperature_at_nodes();

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      map<unsigned int, Real> u;
      map<unsigned int, Real> oldu;
      map<unsigned int, RealGradient> grad_u;

      for (auto&& var : q_var)
      {
        u[var] = 0.0;
        oldu[var] = 0.0;
        grad_u[var] = RealGradient(0);
      }

      // get the solution values at the quadrature point
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (auto&& var : q_var)
        {
          unsigned int ii = id_map[var][i];
          u[var] += phi[i][qp] * Xv.at(var)(ii);
          oldu[var] += phi[i][qp] * oldXv.at(var)(ii);
          grad_u[var] += dphi[i][qp] * Xv.at(var)(ii);
        }
      }

      // add the constant electrochemical potential
      for (auto&& cons : _conservation)
      {
        for (auto&& var : cons.second.carrier_vars)
        {
          // we need this only if the carrier exists in this element
          if (q_var.count(var))
          {
            u[var] += Xv.at(cons.first)(0);
            oldu[var] += oldXv.at(cons.first)(0);
          }
        }
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      double grad_fac = phi0 / x0;
      for (auto&& var : q_var)
      {
        if (var == u_var)
        {
          sc->set_el_potential(phi0 * u[var]);
          sc->set_old_el_potential(phi0 * oldu[var]);
          sc->set_electric_field(-grad_fac * grad_u[var]);
        }
        else
        {
          sc->set_fermi_potential(var, phi0 * u[var]);
          sc->set_old_fermi_potential(var, phi0 * oldu[var]);
          sc->set_grad_fermi(var, grad_fac * grad_u[var]);
        }
      }

      // calculate all local properties
      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_ionized_dopants();
      sc->calculate_mobilities();
      sc->calculate_net_recombination_rates();

      map<unsigned int, long double> dens, R;
      map<unsigned int, double> mu, sigma, tep;

      for (auto&& var : q_var)
      {
        if (var != u_var)
        {
          CarrierProperties* cp = sc->get_carrier_properties(var);
          dens.insert( make_pair(var, sc->get_q_density(var)) );
          R.insert( make_pair(var, sc->get_net_q_recombination_rate(var)) );
          sigma.insert( make_pair(var, cp->get_conductivity() / (mu0 * C0_q)) );
          tep.insert( make_pair(var, cp->get_thermoelectric_power() / phi0) );
        }
      }
      //double Nd = sc->get_ionized_donor_density();
      //double Na = sc->get_ionized_acceptor_density();

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
      for (auto&& var : q_var)
      {

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          // don't need this for constant potentials
          if (dof_indices_var[var].size() > 1)
          {
            for (unsigned int j = 0; j < n_dofs; j++)
            {
              if (var == u_var)
              {
                double laplace_u = J * (dphi[i][qp] * (permittivity * dphi[j][qp]));

                if (coupling & POISSON)
                  Kvv[var].at(var)(i,j) += l2 * laplace_u / scalev.at(var)(i);
              }
              else
              {
                double laplace = J * (dphi[i][qp] * dphi[j][qp]);

                if (coupling & CURRENTS)
                  Kvv[var].at(var)(i,j) += sigma[var] * laplace / scalev.at(var)(i);
              }
            }
          }

          if ((var == u_var) && !(coupling & POISSON))
            Kvv[var].at(var)(i,i) += 1.0;
          if ((var != u_var) && !(coupling & CURRENTS))
          {
            unsigned int ii = id_map[var][i];
            Kvv[var].at(var)(ii,ii) += 1.0;
          }
        }
      }





      //
      // for jacobian compute the other contributions
      //
      if (jacobian != NULL)
      {
        //map<unsigned int, long double> ddens_dphi;
        map<unsigned int, map<unsigned int, long double>> dR;

        for (auto&& vari : q_var)
        {
          if (vari != u_var)
          {
            //long double ddens_dphi_tmp = sc->get_q_density_derivative(vari);
            //ddens_dphi[vari] = ddens_dphi_tmp;

            for (auto&& varj : q_var)
            {
              long double dR_dEf_tmp = sc->get_net_q_recombination_rate_derivatives(vari)[varj];
              long double dR_tmp = dR_dEf_tmp * phi0 / R0;

              //if (R[vari] == 0.0)
              //  dR_tmp = 0.0;

              dR[vari][varj] = dR_tmp;
            }
          }
        }

        for (auto&& var : _conservation)
        {
          for (auto&& vari : var.second.carrier_vars)
          {
            // we need this only if the carrier exists in this element!
            if (q_var.count(vari))
            {
              double ddens = sc->get_q_density_derivative(vari);
              // diagonal part
              Kvv[var.first].at(var.first)(0,0) += J * ddens * phi0 / C0;

              for (unsigned int j = 0; j < n_dofs; j++)
              {
                if (coupling & POISSON)
                  Kvv[var.first].at(u_var)(0,j) -= J * ddens * phi[j][qp] * phi0 / C0;

                unsigned int jj = id_map[vari][j];
                if (coupling & CURRENTS)
                  Kvv[var.first].at(vari)(0,jj) += J * ddens * phi[j][qp] * phi0 / C0;
              }
            }
          }
        }

        map<unsigned int, double> drho;
        drho[u_var] = 0.0;

        for (auto&& var : q_var)
        {
          if (var != u_var)
          {
            drho[var] = sc->get_charge_density_derivative(var) * phi0 / C0;
            drho[u_var] -= drho[var];

            if (params.local_neutrality)
              drho[u_var] = drho[var] = 0.0;
          }
        }

        for (auto&& var : _conservation)
        {
          drho[var.first] = 0.0;
          for (auto&& vari : var.second.carrier_vars)
          {
            // we need this only if the carrier exists in this element!
            if (q_var.count(vari))
              drho[var.first] += sc->get_charge_density_derivative(vari) * phi0 / C0;
          }
        }


        // d(sigma_n)/du * element-jacobian
        // sigma_n is the conductivity of electrons
        // the factor phi_0 comes from the derivative with respect to the potential
        map<unsigned int, double> dsigma;
        map<unsigned int, libMesh::RealGradient> dsigma_grad_u;
        map<unsigned int, libMesh::RealGradient> dsigma_grad_f;

        for (auto&& var : q_var)
        {
          if (var != u_var)
          {
            double der_qFermi;
            libMesh::RealGradient der_grad_u(0);
            libMesh::RealGradient der_grad_f(0);
            sc->get_carrier_properties(var)->
                get_conductivity_and_derivatives(der_qFermi, der_grad_f, der_grad_u);

            double dsigma_tmp = J * phi0 / (mu0 * C0_q) * der_qFermi;
            dsigma.insert( make_pair(var, dsigma_tmp) );

            der_grad_u *= J * phi0 / (mu0 * C0_q) / x0;
            dsigma_grad_u.insert( make_pair(var, der_grad_u) );

            if (dim > 1)
            {
              der_grad_f *= J * phi0 / (mu0 * C0_q) / x0;
              dsigma_grad_f.insert( make_pair(var, der_grad_f) );
            }
          }
        }

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          map<unsigned int, double> lap,
                                    dsigma_x_lap;
          for (auto&& var : q_var)
          {
            if ((var != u_var) && (coupling & CURRENTS) &&
                (dof_indices_var[var].size() > 1))
            {
              lap.insert( make_pair(var, (dphi[i][qp] * grad_u[var]) / scalev.at(var)(i)) );
              dsigma_x_lap.insert ( make_pair(var, dsigma[var]*lap[var]) );
            }
          }

          for (unsigned int j = 0; j < n_dofs; j++)
          {
            // first the dKe_il/dX_j * X_l part
            // (for X_l = u_l we dont get anything, i.e. the
            // contributions to Kuu, Kun, Kup are zero)
            //

            if (_options.exact_newton)
            {

              for (auto&& var : q_var)
              {
                if ((var != u_var) && (coupling & CURRENTS) &&
                    (dof_indices_var[var].size() > 1))
                {
                  double dsigma_x_phi = dsigma_x_lap[var] * phi[j][qp];
                  double dsigma_grad_u_x_dphi = dsigma_grad_u[var] * dphi[j][qp];
                  double dsigma_grad_f_x_dphi = dsigma_grad_f[var] * dphi[j][qp];

                  //if the carrier has zero charge, there is no dependence on the electric potential, hence dsigma_du is zero
                  double q = sc->get_carrier_properties(var)->get_charge();
                  double pot_fac = ( q == 0 ) ? 0.0 : 1.0;

                  if (coupling & POISSON)
                    Kvv.at(var).at(u_var)(i,j) += pot_fac * dsigma_x_phi + dsigma_grad_u_x_dphi * lap[var];

                  Kvv.at(var).at(var)(i,j) += dsigma_grad_f_x_dphi * lap[var] - dsigma_x_phi;



                  double sign = sc->get_carrier_properties(var)->get_charge_sign();
                  // contribution of the Seebeck effect ->residual_derivative
                  double dsigma_e_x_phi_x_P = sign * dsigma[var] * phi[j][qp] * tep[var];
        /*
                  for (unsigned int k = 0; k < n_dofs; k++)
                  {
                    double laplace = dphi[i][qp] * dphi[k][qp];

                    double elem_contrib =
                        dsigma_e_x_phi_x_P * laplace * T_nodes[k] / scalev.at(var)(i);

                    Kvv.at(var).at(var)(i,j) -= elem_contrib;

                    if (coupling & POISSON)
                      Kvv.at(var).at(u_var)(i,j) += elem_contrib;

                  }
         */
                }
              }
            }


            
            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            if (coupling & POISSON)
            {
              Kvv[u_var].at(u_var)(i,j) -= drho[u_var] * phi_i_x_phi_j / scalev.at(u_var)(i);
              //Kvv[u_var].at(u_var)(i,i) -= drho[u_var] * phi_i_x_phi_j / scalev.at(u_var)(i);

              for (auto&& var : q_var)
              {
                if ( (var != u_var) && (coupling & CURRENTS) )
                {
                  unsigned int jj = id_map[var][j];
                  Kvv[u_var].at(var)(i,jj) -= drho[var] * phi_i_x_phi_j / scalev.at(u_var)(i);
                  //Kvv[u_var].at(var)(i,i) -= drho[var] * phi_i_x_phi_j / scalev.at(u_var)(i);
                }
              }
            }

            if (coupling & CURRENTS)
            {
              for (auto&& vari : q_var)
              {
                if (vari != u_var)
                {
                  double sign = sc->get_carrier_properties(vari)->get_charge_sign();
                  unsigned int ii = id_map[vari][i];

                  for (auto&& varj : q_var)
                  {
                    unsigned int jj = id_map[varj][j];
                    Kvv[vari].at(varj)(ii,jj) += sign * dR[vari][varj] * phi_i_x_phi_j / scalev.at(vari)(ii);
                  }

                }
              }
            }
          }

          for (auto&& var : _conservation)
          {
            if (coupling & POISSON)
            {
              Kvv[u_var].at(var.first)(i,0) -=
                  J * drho[var.first] * phi[i][qp] / scalev.at(u_var)(i);

            }

            if (coupling & CURRENTS)
            {
              for (auto&& varj : var.second.carrier_vars)
              {
                // we need this only if the carrier exists in this element!
                if (q_var.count(varj))
                {
                  for (auto&& vari : q_var)
                  {
                    if (vari != u_var)
                    {
                      unsigned int ii = id_map[vari][i];
                      double sign = sc->get_carrier_properties(vari)->get_charge_sign();
                      Kvv[vari].at(var.first)(ii,0) +=
                          J * sign * dR[vari][varj] * phi[i][qp] / scalev.at(vari)(ii);
                    }
                  }
                  //Kvv.at(varj).at(var.first)(i,0) -= dsigma_x_lap[varj];
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

        long double J_x_P0 = J / P0;

        // net recombination rate
        map<unsigned int, long double> J_x_R;
        for (auto&& var : q_var)
        {
          if (var != u_var)
            J_x_R.insert( make_pair(var, J*R[var]/R0) );
        }

        libMesh::RealVectorValue P(sc->get_total_polarization());
        P *= J_x_P0;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (auto&& var : q_var)
          {
            if (var == u_var)
            {
              if (coupling & POISSON)
                Fv.at(var)(i) -= (J_x_rho * phi[i][qp] + (P * dphi[i][qp])) / scalev.at(var)(i);
              else
                Fv.at(var)(i) -= Xv.at(var)(i);
            }
            else
            {
              unsigned int ii = id_map[var][i];
              long double net_recomb = J_x_R[var] * phi[i][qp] / scalev.at(var)(ii);

              double sign = sc->get_carrier_properties(var)->get_charge_sign();

              if (coupling & CURRENTS)
              {
                Fv.at(var)(ii) += sign * net_recomb;

                double sigma_x_P_x_J = J * sign * sigma[var] * tep[var] / scalev.at(var)(ii);
                if (sigma_x_P_x_J != 0)
                {
                  // include Seebeck contribution -> Residual
                  for (unsigned int k = 0; k < n_dofs; k++)
                  {
                    Real laplace = dphi[i][qp] * dphi[k][qp];

                    Fv.at(var)(ii) += sigma_x_P_x_J * laplace * T_nodes[k];
                  }
                }
              }
              else
                Fv.at(var)(ii) -= Xv.at(var)(ii);
            }

          }
        }

        for (auto&& var : _conservation)
        {
          for (auto&& vari : var.second.carrier_vars)
          {
            // we need this only if the carrier exists in this element!
            if (q_var.count(vari))
            {
              double dens = sc->get_q_density(vari);

              Fv.at(var.first)(0) -= J * dens / C0;
            }
          }
        }

      }  //end residual

    } // end loop over quadrature points

    map<unsigned int, vector<double> > dirichlet_jac;
    map<unsigned int, unsigned int > dirichlet_node;
    map<unsigned int, double> dirichlet_res;

    // now loop over the element sides to find boundary elements
    // and to include von Neumann and mixed type boundary conditions
    //
    // NOTE 1:
    // we dont apply BC for nabla(Ef) but for the particle
    // flux sigma * nabla(Ef)
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
            u.insert( make_pair(var, 0) );
            grad_u.insert( make_pair(var, libMesh::RealGradient(0)) );
          }

          // get the solution values at the quadrature point
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            for (auto var : q_var)
            {
              unsigned int ii = id_map[var][i];
              u[var] += phi_face[i][qp] * Xv.at(var)(ii);
              grad_u[var] += dphi_face[i][qp] * Xv.at(var)(ii);
            }
          }

          // add the constant electrochemical potential
          for (auto&& cons : _conservation)
          {
            u[cons.second.carrier_vars[0]] += Xv.at(cons.first)(0);
            oldu[cons.second.carrier_vars[0]] += oldXv.at(cons.first)(0);
          }

          // prepare for calculating local properties
          sc->set_coordinates(q_point_face[qp]);

          double grad_fac = phi0 / x0;
          for (auto var : q_var)
          {
            //cout<<"u"<<var<<" = "<<u[var]<<endl;
            if (var == u_var)
            {
              sc->set_el_potential(phi0 * u[var]);
              sc->set_electric_field(-grad_fac * grad_u[var]);
            }
            else
            {
              sc->set_fermi_potential(var, phi0 * u[var]);
              sc->set_grad_fermi(var, grad_fac * grad_u[var]);
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
            if (var == u_var)
            {
              sm->set_el_potential(phi0 * u[var]);
              sm->set_electric_field(-grad_fac * grad_u[var]);
            }
            else
            {
              sm->set_fermi_potential(var, phi0 * u[var]);
              sm->set_grad_fermi(var, grad_fac * grad_u[var]);
            }
          }

          sm->compute();

          // the jacobian x weight x scaling
          double J = JxW_face[qp];


          // contribution to the jacobian
          if (jacobian != NULL)
          {
            // for Dirichlet DOFs we do not add anything
            vector<double> scale_v(_carriers.size() + 1, 0.0);
            vector<vector<double>> deriv_v(_carriers.size() + 1, vector<double>(_carriers.size() + 1, 0));
            for (auto&& var : q_var)
            {
              scale_v[var] = (var==u_var) ?  J * phi0 / x0 / C0 : J * phi0 / (x0 * R0);
              if (sm->get_type(var) == DDInterfaceModel::DIRICHLET)
                scale_v[var] = 0;

              vector<double> deriv_tmp = sm->get_jacobian_row(var);
              deriv_v[var] = deriv_tmp;
            }


            for (unsigned int i = 0; i < n_dofs; i++)
            {
              vector<double> fac(_carriers.size() + 1, 0.0);

              if (elem->is_node_on_side(i, s))
              {
                unsigned int n_var = 0;
                for (auto&& var : q_var)
                {
                  unsigned int ii = id_map[var][i];
                  fac[var] = scale_v[var] / scalev.at(var)(ii);

                  if (sm->get_type(var) == DDInterfaceModel::DIRICHLET)
                  {
                    dirichlet_jac[n_var + ii] = deriv_v[var];
                    dirichlet_node[n_var + ii] = i;
                  }

                  n_var += dof_indices_var[var].size();
                }
              }

              for (unsigned int j = 0; j < n_dofs; j++)
              {
                Real phi_i_x_phi_j = phi_face[i][qp] * phi_face[j][qp];

                for (auto&& vari : q_var)
                {
                  if ((vari == u_var) && (coupling & POISSON))
                  {
                    for (auto&& varj : q_var)
                    {
                      unsigned int jj = id_map[varj][j];
                      if ((varj == u_var) || ((varj != u_var) && (coupling & CURRENTS)))
                        Kvv.at(vari).at(varj)(i,jj) -= fac[vari] * deriv_v[vari][varj] * phi_i_x_phi_j;
                    }
                  }
                  else if ((vari != u_var) && (coupling & CURRENTS))
                  {
                    unsigned int ii = id_map[vari][i];

                    const char ct = sc->get_carrier_properties(vari)->get_carrier_type();
                    double sign = (ct == 'e') ? 1.0 : -1.0;
                    for (auto&& varj : q_var)
                    {
                      unsigned int jj = id_map[varj][j];

                      if ( ((varj == u_var) && (coupling & POISSON)) || (varj != u_var) )
                        Kvv.at(vari).at(varj)(ii,jj) -= sign * fac[vari] * deriv_v[vari][varj] * phi_i_x_phi_j;
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
            vector<double> value_v(_carriers.size() + 1, 0.0);

            const vector<double>& coeff_a = sm->get_a();
            const vector<double>& coeff_g = sm->get_g();

            for (auto&& var : q_var)
            {
              value_v[var] = (var == u_var) ? J / (x0 * C0) : J / (x0 * R0);
              if (sm->get_type(var) != DDInterfaceModel::DIRICHLET)
                value_v[var] *= (coeff_g[var] - coeff_a[var] * u[var] * phi0);
              else
                value_v[var] = coeff_g[var] / phi0;

              //cout<<"g"<<var<<" = "<<coeff_g[var]<<" a"<<var<<" = "<<coeff_a[var]<<" v"<<var<<" = "<<u[var]<<endl;
            }



            for (unsigned int i = 0; i < n_dofs; i++)
            {
              unsigned int n_var = 0;
              for (auto&& var : q_var)
              {
                unsigned int ii = id_map[var][i];
                double scale_v = 1.0 / scalev.at(var)(ii);

                if (elem->is_node_on_side(i, s))
                {
                  if (sm->get_type(var) == DDInterfaceModel::DIRICHLET)
                  {
                    dirichlet_res[n_var + ii] = value_v[var];
                    scale_v = 0;
                  }
                }
                n_var += dof_indices_var[var].size();

                if ((var == u_var) && (coupling & POISSON))
                  Fv.at(var)(i) -= value_v[var] * phi_face[i][qp] * scale_v;
                else if ((var != u_var) && (coupling & CURRENTS))
                {
                  const char ct = sc->get_carrier_properties(var)->get_carrier_type();
                  double sign = (ct == 'e') ? 1.0 : -1.0;
                  Fv.at(var)(ii) -= sign * value_v[var] * phi_face[i][qp] * scale_v;
                }

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
                Fv.at(u_var)(i) -= value_u * phi_face[i][qp] / scalev.at(u_var)(i);
              }
            }

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


        unsigned int j = dirichlet_node[i];
        unsigned int ctr = 0;
        for (auto&& var : q_var)
        {
          unsigned int jj = id_map[var][j];
          Ke(i,ctr + jj) = -(it->second)[var];
          ctr += dof_indices_var[var].size();
        }
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
      i
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
    //jacobian->print_matlab("J.m");
    //exit(0);

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
      jacobian->print_matlab("J.m");
      //ostringstream os;
      //os << "_" << __private_counter << ".m";
      //ostringstream ms;
      //ms << "writing " << "J" << os.str() << "\n";
      //Messages::info(ms.str());
      //jacobian->print_matlab("J" + os.str());
      //sysmat.print_matlab("precond" + os.str());
    }
    if (coupling & ELECTRONS) __private_counter++;
    //if (__private_counter == 2) exit(0);
    */
  }
  else
  {
    // the following is done only on node 0, otherwise we would
    // add up the constant term many times
    //if (this->get_communicator().rank() == 0)
    //{
      for (auto&& var : _conservation)
      {
        vector<dof_id_type> scalars;
        //dof_map.SCALAR_dof_indices(scalars, var.first);
        dof_map.local_variable_indices(scalars, mesh, var.first);

        double scale = x0 * C0;
        switch (dim)
        {
          case 3:
            scale *= x0;

          case 2:
            scale *= x0;

          default:
            break;
        }
        if (!scalars.empty())
          residual->add(scalars[0], var.second.conserved_number / scale);
      }
    //}

    residual->close();
    //residual->print_matlab("r.m");
    //exit(0);

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

  vector<double> results(3 * nn, 0.0);

  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("el_qw");
  unsigned int ep_var = system->variable_number("hl_qw");

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
      if (dof_indices_en.size() > 0)
        results[id + 1] = vec(dof_indices_en[n]);
      if (dof_indices_ep.size() > 0)
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

