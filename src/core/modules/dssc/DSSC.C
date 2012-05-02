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
#include "TiberLinearSystem.h"
#include "SolveFailedException.h"
#include "Messages.h"


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


// C++ includes

using namespace std;


//
// Module interface
//

TIBER_MODULE(DSSC, dssc)




DSSC*
DSSC::_this;




DSSC::DSSC(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true),
    _poisson_only(false),
    _scaling_type(Scaling::UNITS)
{
}




DSSC::~DSSC(void)
{
  cleanup_solver();
}




PhysicalModel*
DSSC::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
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

  // we calculate in cm!
  double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  unsigned int dim = get_mesh().mesh_dimension();

  _cond_scaling.n = 1;
  _cond_scaling.I = 1;
  _cond_scaling.I3 = 1;
  _cond_scaling.C = 1;

  /*
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
  */

  // the scaling parameters should never be zero
  // they are in any case positive, so it will
  // always find the maximum
  double x0 = -1;
  double phi0 = SimulationOptions::T * Constants::k_B;
  double mu0 = 1;
  double C0 = 1;
  double eps0 = -1;

  _cation_amount = 0.0;
  _iodine_amount = 0.0;

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

    sc->reinit(elem);
    sc->set_coordinates(elem->centroid());
    sc->set_potentials(0.0);
    sc->set_electric_field(RealGradient(0));
    sc->set_grad_fermi_n(RealGradient(0));
    sc->set_grad_fermi_I(RealGradient(0));
    sc->set_grad_fermi_I3(RealGradient(0));
    sc->set_grad_fermi_C(RealGradient(0));

    sc->calculate_densities();
    sc->calculate_traps();
    sc->calculate_equilibrium_traps();
    sc->calculate_net_recombination_rate();

    // element volume in cm
    double volume = elem->volume();
    switch (dim)
    {
      case 3:
        volume *= mesh_units;
      case 2:
        volume *= mesh_units;
      default:
        volume *= mesh_units;
    };
    _cation_amount += sc->get_equilibrium_concentrations().C * volume;
    _iodine_amount += sc->get_equilibrium_concentrations().I3 * volume;
    _iodine_amount += sc->get_equilibrium_concentrations().I * volume / 3.0;

    double mu = sc->get_mobility_n();
    double C = sc->get_equilibrium_concentrations().n;
    //C = (C > C1) ? C : C1;
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.n = (mu * C > _cond_scaling.n) ? mu * C : _cond_scaling.n;

    mu = sc->get_mobility_I();
    C = sc->get_equilibrium_concentrations().I;
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.I = (mu * C > _cond_scaling.I) ? mu * C : _cond_scaling.I;

    mu = sc->get_mobility_I3();
    C = sc->get_equilibrium_concentrations().I3;
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.I3 = (mu * C > _cond_scaling.I3) ? mu * C : _cond_scaling.I3;

    mu = sc->get_mobility_C();
    C = sc->get_equilibrium_concentrations().C;
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.C = (mu * C > _cond_scaling.C) ? mu * C : _cond_scaling.C;


    double eps = sc->get_relative_permittivity();
    eps0 = (eps0 > eps) ? eps0 : eps;
  }

  switch (type)
  {
    default: // UNITS

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

  //mu0 = 1;
  //C0 = 1;
  //eps0 = 1;
  x0 = 1;
  mesh_units = 1;
  phi0 = 1;
  get_scaling().set_potential_scaling(phi0);
  get_scaling().set_length_scaling(x0 * mesh_units);
  get_scaling().set_mobility_scaling(mu0);
  get_scaling().set_density_scaling(C0);

  ostringstream os;
  os << "total amount cation: " << _cation_amount << endl;
  os << "total amount iodine: " << _iodine_amount << endl;
  Messages::info(os.str());
}




void
DSSC::compute_scaling_only(Scaling::ScalingType type)
{

  // we calculate in cm!
  //double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  //get_scaling().set_calc_mesh_units(mesh_units);

  const MeshBase& mesh = get_mesh();

  _cond_scaling.n = 1;
  _cond_scaling.I = 1;
  _cond_scaling.I3 = 1;
  _cond_scaling.C = 1;

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

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const NumericVector<Number>& solution = get_solution_vector();


  const DofMap& dof_map = system->get_dof_map();

    // the scaling parameters to scale back the result
    phi0 = get_scaling().get_potential_scaling();
    C0 = get_scaling().get_density_scaling();

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
    const unsigned int dim = mesh.mesh_dimension();

    FEType fe_type = system->variable_type(u_var);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    QGauss qrule(dim, libMeshEnums::CONSTANT);
    fe->attach_quadrature_rule(&qrule);



  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();


    fe->reinit(elem);
    const vector<vector<Real> >& phi = fe->get_phi();
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();

      dof_map.dof_indices(elem, dof_indices_u, u_var);
      dof_map.dof_indices(elem, dof_indices_en, en_var);
      dof_map.dof_indices(elem, dof_indices_eI, eI_var);
      dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
      dof_map.dof_indices(elem, dof_indices_eC, eC_var);

      assert(_device->get_material(elem->subdomain_id()) != NULL);
      DSSCModel* sc =
        dynamic_cast<DSSCModel*>(
            _device->get_material(elem->subdomain_id())->get_model(get_id()));

      sc->reinit(elem);
      // prepare for calculating local properties
      sc->set_coordinates(elem->centroid());

      unsigned int n_dofs = dof_indices_u.size();
      // get the solution values at the centroid
      double u  = 0.0;
      double en = 0.0;
      double eI = 0.0;
      double eI3 = 0.0;
      double eC = 0.0;
      RealGradient e_field(0);
      RealGradient grad_en(0);
      RealGradient grad_eI(0);
      RealGradient grad_eI3(0);
      RealGradient grad_eC(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][0] * solution(dof_indices_u[i]);
        en += phi[i][0] * solution(dof_indices_en[i]);
        eI += phi[i][0] * solution(dof_indices_eI[i]);
        eI3 += phi[i][0] * solution(dof_indices_eI3[i]);
        eC += phi[i][0] * solution(dof_indices_eC[i]);

        e_field += dphi[i][0] * solution(dof_indices_u[i]);
        grad_en += dphi[i][0] * solution(dof_indices_en[i]);
        grad_eI += dphi[i][0] * solution(dof_indices_eI[i]);
        grad_eI3 += dphi[i][0] * solution(dof_indices_eI3[i]);
        grad_eC += dphi[i][0] * solution(dof_indices_eC[i]);
      }
      e_field *= -phi0;
      grad_en *= phi0;
      grad_eI *= phi0;
      grad_eI3 *= phi0;
      grad_eC *= phi0;


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);

      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_equilibrium_traps();
      sc->calculate_net_recombination_rate();

    double mu = sc->get_mobility_n();
    double C = sc->get_density_n();
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.n = (mu * C > _cond_scaling.n) ? mu * C : _cond_scaling.n;

    mu = sc->get_mobility_I();
    C = sc->get_density_I();
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.I = (mu * C > _cond_scaling.I) ? mu * C : _cond_scaling.I;

    mu = sc->get_mobility_I3();
    C = sc->get_density_I3();
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.I3 = (mu * C > _cond_scaling.I3) ? mu * C : _cond_scaling.I3;

    mu = sc->get_mobility_C();
    C = sc->get_density_C();
    mu0 = (mu0 > mu) ? mu0 : mu;
    C0 = (C0 > C) ? C0 : C;
    _cond_scaling.C = (mu * C > _cond_scaling.C) ? mu * C : _cond_scaling.C;


    double eps = sc->get_relative_permittivity();
    eps0 = (eps0 > eps) ? eps0 : eps;
  }

  switch (type)
  {
    default: // UNITS

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

  x0 = 1;
  double mesh_units = 1;
  phi0 = 1;
  get_scaling().set_potential_scaling(phi0);
  get_scaling().set_length_scaling(x0 * mesh_units);
  get_scaling().set_mobility_scaling(mu0);
  get_scaling().set_density_scaling(C0);

}



void
DSSC::get_OC_values(void)
{
/*
  // references for nicer code
  const MeshBase& mesh = get_mesh();
  EquationSystems& eq_sys = get_equation_systems();
  TiberNonlinearSystem& system = static_cast<TiberNonlinearSystem&>(
      eq_sys.get_system(get_equation_system_name()));

  const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  const NumericVector<Number>& solution = system.get_solution_vector();

  const unsigned int dim = mesh.mesh_dimension();

  const Scaling& scaling = get_scaling();
  const double x0 = scaling.get_length_scaling();
  const double phi0 = scaling.get_potential_scaling();

  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_n");
  const unsigned int eI_var = system.variable_number("fermi_I");
  const unsigned int eI3_var = system.variable_number("fermi_I3");
  const unsigned int eC_var = system.variable_number("fermi_C");

  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;



  MeshBase::const_element_iterator el =
                                  mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_local_elements_end();

  // loop over all active elements
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = elem->top_parent();

    ID subdomain = elem->subdomain_id();

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(top_parent, s);

      // is this a boundary?
      if (environment.is_boundary(side))
      {
        Boundary* boundary = environment.get_boundary(side);

        DSSCContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<DSSCContact*>(
              boundary->get_boundary_properties(get_id()));

            if ( (!contact->is_cathode()) && (!contact->is_gate()) )
        if ((contact != NULL) && (!contact->is_cathode() && !contact->is_gate()))
        {

          dof_map.dof_indices(elem, dof_indices);
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

          double u  = phi0 * solution(dof_indices_u[s]);
          double en = phi0 * solution(dof_indices_en[s]);
          double eI = phi0 * solution(dof_indices_eI[s]);
          double eI3 = phi0 * solution(dof_indices_eI3[s]);
          double eC = phi0 * solution(dof_indices_eC[s]);


          // prepare for calculating local properties
          sc->set_coordinates(elem->point(s));

          sc->set_potentials(u, en, eI, eI3, eC);

          // calculate all local properties
          sc->calculate_densities();

          //double n_e = sc->get_density_n();
//          double n_I = sc->get_density_I();
          double n_I = sc->get_density_e();
          double nIdark = sc->get_equilibrium_concentrations().I;
          double n_I3 = sc->get_density_I3();
          double nI3dark = sc->get_equilibrium_concentrations().I3;
          //double n_C = sc->get_density_C();

          contact->set_values(n_I, nIdark, n_I3, nI3dark);

        }
      }
    }
  }

*/
}






void
DSSC::reset_solver(void)
{
  if (!_rebuild_eq_system)
  {
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
    //get_OC_values();
    compute_scaling_only(_scaling_type);
    //build_local_scaling();
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

  Messages::newline();
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
    Messages::info(os.str());
  }

  for ( ; it != end; ++it)
  {

    ostringstream os;
    os << setprecision(6);
    DSSCContact* cnt =
      static_cast<DSSCContact*>(it->first->get_boundary_properties(get_id()));
    os << it->first->get_name();
    os.width(width - os.tellp());
    os << "";
    os << (cnt->get_potential()/it->first->get_area_factor()) * it->second;
    os.width(2 * width - os.tellp());
    os << "";
    os << it->second;
    Messages::info(os.str());

  }

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

}



void
DSSC::guess_equilibrium(void)
{

}



void
DSSC::do_print_info(void)
{
  parse_const_options();
  parse_options();

  if (do_local_scaling_)
    Messages::info("using local scaling");


}





void
DSSC::parse_const_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();

  //do_local_scaling_ = opts.get_option("local_scaling", true);
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




  ModelOptions& solveropts = get_solver_options();
  if (solveropts.get_option("absolute_tolerance", -1.0) < 0)
    solveropts["absolute_tolerance"] = "1e-15";



  // the coupled DD system
  create_equation_system("nonlinear");
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();
  //TiberNonlinearSystem& system =
  //  *TiberNonlinearSystem::create(equation_systems,
  //    get_equation_system_name(), get_solver_options());

  system.attach_assembly_routine(assemble_system);


  system.add_variable("potential", libMeshEnums::FIRST);
  system.add_variable("fermi_n", libMeshEnums::FIRST);
  system.add_variable("fermi_I", libMeshEnums::FIRST);
  system.add_variable("fermi_I3", libMeshEnums::FIRST);
  system.add_variable("fermi_C", libMeshEnums::FIRST);

  // finally initialize the newly created system
  system.init();

/*  if (EIS .eq. true)
  {
    // the coupled DD system
    create_equation_system("linear");
    TiberLinearSystem& system_frequency = get_equation_system<TiberLinearSystem>(1);

    system_frequency.attach_assembly_routine(assemble_system);

    system_frequency.add_variable("potential_R", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_n_R", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_I_R", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_I3_R", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_C_R", libMeshEnums::FIRST);
    system_frequency.add_variable("potential_I", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_n_I", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_I_I", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_I3_I", libMeshEnums::FIRST);
    system_frequency.add_variable("fermi_C_I", libMeshEnums::FIRST);

    system_frequency.init();
    }
*/
  _rebuild_eq_system = false;

}





void
DSSC::do_init(void)
{

  _device = &get_environment().get_device();

  find_dirichlet_nodes();
  find_internal_boundary_nodes();

  parse_const_options();

  rebuild_equation_system();

  // prepare the _boundary_currents
  // we will rely on the fact that it contains an entry for every boundary
  SimulationEnvironment::BoundaryIterator
    it(get_environment().boundaries_begin());
  const SimulationEnvironment::BoundaryIterator
    end(get_environment().boundaries_end());
  for ( ; it != end; ++it)
  {
      BoundaryProperties* bd = (*it)->get_boundary_properties(get_id());
      if (bd != NULL)
      {
        DSSCContact* contact = dynamic_cast<DSSCContact*>(bd);
        if (contact != NULL)
        {
          _boundary_currents[(*it)] = 0.0;
          _voltages[(*it)] = 0.0;
        }
      }
    }

}





NumericVector<double>&
DSSC::solution_vector(void)
{
  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  return system.get_solution_vector();
}




void
DSSC::do_newton(void)
{

  EquationSystems& es = get_equation_systems();

  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();


  system.set_options(get_solver_options());
  system.solve();
}





void
DSSC::find_dirichlet_nodes(void)
{
  SimulationEnvironment& env = get_environment();

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
        Boundary* bd = env.get_boundary(side);
        if (bd == NULL)
          continue;


        DSSCContact* contact = NULL;
        contact = dynamic_cast<DSSCContact*>(
            bd->get_boundary_properties(get_id()));

        if (contact != NULL)
        {
          AutoPtr<Elem> side_el(elem->build_side(s));
          for (unsigned int i = 0; i < side_el->n_nodes(); i++)
            _dirichlet_nodes[side_el->get_node(i)] = bd;
        }
      }
    }
  }

}




void
DSSC::calculate_currents_rstf(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;

  // reset currents
  ContactData::iterator it =
    _boundary_currents.begin();
  for ( ; it != _boundary_currents.end(); ++it)
    (*it).second = 0.0;

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

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
  const unsigned int en_var = system->variable_number("fermi_n");
  const unsigned int eI_var = system->variable_number("fermi_I");
  const unsigned int eI3_var = system->variable_number("fermi_I3");
  const unsigned int eC_var = system->variable_number("fermi_C");

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
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;


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
    dof_map.dof_indices(elem, dof_indices_eI, eI_var);
    dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
    dof_map.dof_indices(elem, dof_indices_eC, eC_var);

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
      Real eI = 0.0;
      Real eI3 = 0.0;
      Real eC = 0.0;
      RealGradient dEfn(0);
      RealGradient dEfI(0);
      RealGradient dEfI3(0);
      RealGradient dEfC(0);
      RealGradient e_field(0);
      RealGradient dT(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        eI += phi[i][qp] * solution(dof_indices_eI[i]);
        eI3 += phi[i][qp] * solution(dof_indices_eI3[i]);
        eC += phi[i][qp] * solution(dof_indices_eC[i]);

        dEfn += dphi[i][qp] * solution(dof_indices_en[i]);
        dEfI += dphi[i][qp] * solution(dof_indices_eI[i]);
        dEfI3 += dphi[i][qp] * solution(dof_indices_eI3[i]);
        dEfC += dphi[i][qp] * solution(dof_indices_eC[i]);

        //dT += dphi[i][qp] * T_nodes[i];

        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }

      // prepare for calculating local properties
      sc->set_coordinates(elem->centroid());


      sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);

      sc->set_electric_field(e_field);
      sc->set_grad_fermi_n(dEfn);
      sc->set_grad_fermi_I(dEfI);
      sc->set_grad_fermi_I3(dEfI3);
      sc->set_grad_fermi_C(dEfC);

      sc->calculate_densities();
      sc->calculate_traps();
      sc->calculate_equilibrium_traps();
      sc->calculate_net_recombination_rate();

      // we put the minus here for convenience
      double sigma_I = -Constants::e * sc->get_mobility_I() * sc->get_density_I() ;
      double sigma_I3 = -Constants::e * sc->get_mobility_I3() * sc->get_density_I3() ;
      double sigma_n = -Constants::e * sc->get_mobility_n() * sc->get_density_n() ;
      double ne = sc->get_density_n() ;
      //double sigma_n = -Constants::e * sc->get_mobility_n() * pow(ne,1.4);

      //RealGradient j(JxW[qp] * phi0 *
      //    (sigma_n * (dEfn + Pn * dT) + sigma_h * (dEfp + Pp * dT)));
      RealGradient j(JxW[qp] * phi0 *
          (sigma_n * dEfn + sigma_I * dEfI + sigma_I3 * dEfI3));

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        Boundary* boundary = node_ids[n];
        if (boundary != NULL)
        {
          DSSCContact* contact = dynamic_cast<DSSCContact*>(
                boundary->get_boundary_properties(get_id()));
          //if (contact->is_real_contact())
            _boundary_currents[boundary] += j * dphi[n][qp];
        }

      }
    } // end loop over quadrature points
  } // end loop over elements

}





void
DSSC::build_local_scaling(void)
{

 /* TiberNonlinearSystem* system =
    &get_equation_systems().get_system<TiberNonlinearSystem>(
        get_equation_system_name());

  const NumericVector<Number>& solution = system->get_solution_vector();

  // aliases for nicer code
  const Device& device = *(_device);
  const MeshBase& mesh = get_mesh();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int dim = mesh.mesh_dimension();

  //const Options& params = get_options();

  // the scaling parameters to scale back the result
  const Scaling& scaling = get_scaling();
  const double phi0 = scaling.get_potential_scaling();
  const double x0 = scaling.get_length_scaling();
  const double C0 = scaling.get_density_scaling();
  const double mu0 = scaling.get_mobility_scaling();
  const double l2 = scaling.get_lambda_squared() * Constants::e0 * 1e-2;


  // numeric ids corresponding to the variables
  const unsigned int u_var = system->variable_number("potential");
  const unsigned int en_var = system->variable_number("fermi_n");
  const unsigned int eI_var = system->variable_number("fermi_I");
  const unsigned int eI3_var = system->variable_number("fermi_I3");
  const unsigned int eC_var = system->variable_number("fermi_C");

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
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;


  local_scaling_.clear();
  {
    MeshBase::const_element_iterator it =
      mesh.active_elements_begin();
    const MeshBase::const_element_iterator end =
      mesh.active_elements_end();

    for ( ; it != end; ++it)
      for (unsigned int n = 0; n < (*it)->n_nodes(); n++)
        local_scaling_[(*it)->get_node(n)] = vector<double>(5, 0.0);
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
    dof_map.dof_indices(elem, dof_indices_eI, eI_var);
    dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
    dof_map.dof_indices(elem, dof_indices_eC, eC_var);

    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));


    fe->reinit(elem);

    sc->reinit(elem);

    assert(elem->n_nodes() == dof_indices_u.size());

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      unsigned int n_dofs = dof_indices_u.size();
      // get the solution values at the centroid
      Real u  = 0.0;
      Real en = 0.0;
      Real eI = 0.0;
      Real eI3 = 0.0;
      Real eC = 0.0;
      RealGradient dEfn(0);
      RealGradient dEfI(0);
      RealGradient dEfI3(0);
      RealGradient dEfC(0);
      RealGradient e_field(0);
      RealGradient dT(0);
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u  += phi[i][qp] * solution(dof_indices_u[i]);
        en += phi[i][qp] * solution(dof_indices_en[i]);
        eI += phi[i][qp] * solution(dof_indices_eI[i]);
        eI3 += phi[i][qp] * solution(dof_indices_eI3[i]);
        eC += phi[i][qp] * solution(dof_indices_eC[i]);

        dEfn += dphi[i][qp] * solution(dof_indices_en[i]);
        dEfI += dphi[i][qp] * solution(dof_indices_eI[i]);
        dEfI3 += dphi[i][qp] * solution(dof_indices_eI3[i]);
        dEfC += dphi[i][qp] * solution(dof_indices_eC[i]);
        e_field += dphi[i][qp] * solution(dof_indices_u[i]);
      }

      // prepare for calculating local properties
      sc->set_coordinates(q_point[qp]);

      sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);

      sc->set_electric_field(e_field);
      sc->set_grad_fermi_n(dEfn);
      sc->set_grad_fermi_I(dEfI);
      sc->set_grad_fermi_I3(dEfI3);
      sc->set_grad_fermi_C(dEfC);

      sc->calculate_densities();
      sc->calculate_net_recombination_rate();

      // the jacobian x weight x scaling
      double J = JxW[qp];

      double mu_I = sc->get_mobility_I();
      double mu_I3 = sc->get_mobility_I3();
      double mu_C = sc->get_mobility_C();
      double mu_n = sc->get_mobility_n();

      double sigma_I = mu_I * sc->get_density_I()* J;
      double sigma_I3 = mu_I3 * sc->get_density_I3() * J;
      double sigma_C = mu_C * sc->get_density_C() * J;
      double sigma_n = mu_n * sc->get_density_n() * J;

      double dn_dphi = sc->get_density_derivative_n();
      double dI_dphi = sc->get_density_derivative_I();
      double dI3_dphi = sc->get_density_derivative_I3();
      double dC_dphi = sc->get_density_derivative_C();
      double epsilon = sc->get_relative_permittivity();
      double l2_eps = J * l2 * epsilon;

      double R = 0.0;
      R = sc->get_net_recombination_rate();

      double drho = dC_dphi - (dn_dphi + dI_dphi + dI3_dphi);
        drho *= phi0;

      double dR[4];

      if (R == 0.0)
        dR[0] = dR[1] = dR[2] = dR[3] = 0.0;
      else
      {
        double dR_dn = sc->get_net_recombination_rate_derivatives()[0];
        double dR_dI = sc->get_net_recombination_rate_derivatives()[1];
        double dR_dI3 = sc->get_net_recombination_rate_derivatives()[2];

        dR[1] = -dR_dn * dn_dphi * phi0 * J;
        dR[2] = -dR_dI * dI_dphi * phi0 * J;
        dR[3] = -dR_dI3 * dI3_dphi * phi0 * J;
        dR[0] = -(dR[1] + dR[2] + dR[3]);
      }


      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //local_scaling_[elem->get_node(i)][0] += sigma_e * phi[i][qp];
        //local_scaling_[elem->get_node(i)][1] += sigma_h * phi[i][qp];
        local_scaling_[elem->get_node(i)][0] += sigma_n * (dphi[i][qp] * dphi[i][qp])
                                                  - dR[1] * (phi[i][qp] * phi[i][qp]);
        local_scaling_[elem->get_node(i)][1] += sigma_I * (dphi[i][qp] * dphi[i][qp])
                                                  + 1.5 * dR[2] * (phi[i][qp] * phi[i][qp]);
        local_scaling_[elem->get_node(i)][2] += sigma_I3 * (dphi[i][qp] * dphi[i][qp])
                                                  - 0.5 * dR[3] * (phi[i][qp] * phi[i][qp]);
        local_scaling_[elem->get_node(i)][3] += sigma_C * (dphi[i][qp] * dphi[i][qp]);
        local_scaling_[elem->get_node(i)][4] += l2_eps * (dphi[i][qp] * dphi[i][qp])
                                                  - J * drho * (phi[i][qp] * phi[i][qp]);
      }


    } // end loop over quadrature points
  } // end loop over elements
*/
}



NumericVector<double>&
DSSC::do_get_solution_vector(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  return system.get_solution_vector();
}



void
DSSC::calculate_currents(void)
{
  calculate_currents_rstf();
}



double
DSSC::do_maximum_norm_of_difference(ID id)
{
  double norm = SimulationInterface::do_maximum_norm_of_difference(id);

  return norm * get_scaling().get_potential_scaling();
}



void
DSSC::find_internal_boundary_nodes(void)
{
  MeshBase& mesh = get_mesh();
  MeshBase::element_iterator it = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for ( ; it != end; ++it)
  {
    const Elem* el = *it;

    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          get_physical_model(el->subdomain_id()));

    // we are only interested in boundaries between semiconductor/dielectric
    if (sc->is_TiO2())
    {
      for (unsigned s = 0; s < el->n_sides(); s++)
      {
        if (get_environment().is_inner_boundary(ElementSide(el, s)))
        {
          // get the model of the neighbor element
          DSSCModel* scn =
            dynamic_cast<DSSCModel*>(
                get_physical_model(el->neighbor(s)->subdomain_id()));


          // if neighbor is not dielectric we record it
          if (!scn->is_TiO2())
          {
            AutoPtr<Elem> side(el->build_side(s));
            for (unsigned int i = 0; i < side->n_nodes(); i++)
              _internal_boundary_nodes.insert(side->get_node(i));
          }
        }
      }
    }
  }
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
  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const unsigned int dim = mesh.mesh_dimension();

  const Device& device = *_device;
  const SimulationEnvironment& environment = get_environment();

  double pot = 0.0;

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
  double C0_e = _cond_scaling.n;
  //double C0_e = _cond_scaling.C;
  // scaling for I
  double C0_I = _cond_scaling.I;
  // scaling for I3
  double C0_I3 = _cond_scaling.I3;
  // scaling for C
  double C0_C = _cond_scaling.C;
  // scaling for I + I3
  double C0_tot = C0;
  // scaling C
  double scaling_C = 1.0;
  // scaling I
  double scaling_I = 1.0;
  // scaling I3
  double scaling_I3 = 1.0;
  // scaling all charges
  double scaling_tot = 1.0;


  if (do_local_scaling_)
  {
    scaling_C = _cond_scaling.C;
    scaling_I = _cond_scaling.I;
    scaling_I3 = _cond_scaling.I3;
    scaling_tot = scaling.get_density_scaling();
    C0_tot = C0_e = C0_I = C0_I3 = C0_C = 1.0;
  }



  // scaling for recombination rates
  double tmp = phi0 / (x0 * x0);
  double R0_e = C0_e * tmp;
  double R0_I = C0_I * tmp ;
  double R0_I3 = C0_I3 * tmp;


  const DofMap& dof_map = system.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var = system.variable_number("potential");
  const unsigned int en_var = system.variable_number("fermi_n");
  const unsigned int eI_var = system.variable_number("fermi_I");
  const unsigned int eI3_var = system.variable_number("fermi_I3");
  const unsigned int eC_var = system.variable_number("fermi_C");

  FEType fe_type = system.variable_type(u_var);

  libMeshEnums::Order integration_order = libMeshEnums::FIFTH;
  //libMeshEnums::Order integration_order = libMeshEnums::FIRST;

  // the finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, integration_order);
  //QTrap qrule(dim);
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


  const Node* n_cat = NULL;
  double tot_cat = 0;
  AutoPtr<NumericVector<Number> > cons_cat = x.clone();
  cons_cat->zero();

  const Node* n_iodine = NULL;
  double tot_iodine = 0;
  AutoPtr<NumericVector<Number> > cons_iodine = x.clone();
  cons_iodine->zero();

  set<unsigned int> inner_boundary_nodes;

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

    vector<vector<double> > local_scaling(elem->n_nodes(), vector<double>(5, 1));
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
      sc->calculate_traps();
      sc->calculate_equilibrium_traps();
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
      double sigma_n = (mu_n * n_e ) / C0_e;
      double sigma_I = (mu_I * n_I ) / C0_I;
      double sigma_I3 = (mu_I3 * n_I3 ) / C0_I3;
      double sigma_C = (mu_C * n_C ) / C0_C;

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
          double laplace = J * (dphi[i][qp] * dphi[j][qp]);

          Kuu(i,j) += l2_eps * laplace / local_scaling[i][4];

          if (!poisson_only())
          {
            if (is_TiO2)
              Knn(i,j) += sigma_n * laplace / local_scaling[i][0];

            if (is_electrolyte)
            {
              KII(i,j) += sigma_I * laplace / local_scaling[i][1];

              KI3I3(i,j) += sigma_I3 * laplace / local_scaling[i][2];

              KCC(i,j) += sigma_C * laplace / local_scaling[i][3];
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
        double dn_dphi_exp_DOS = sc->get_density_derivative_n_exp_DOS();

        double drho = dC_dphi - (dn_dphi_exp_DOS + dI_dphi + dI3_dphi);
        drho *= phi0 / C0_tot;


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
        double dsigma_n = J * phi0 / C0_e * mu_n * dn_dphi;
        double dsigma_I = J * phi0 / C0_I * mu_I * dI_dphi;
        double dsigma_I3 = J * phi0 / C0_I3 * mu_I3 * dI3_dphi;
        double dsigma_C = J * phi0 / C0_C * mu_C * dC_dphi;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          //double lap_e = (dphi[i][qp] * grad_en) / local_scaling[i][0];
          //double lap_I = (dphi[i][qp] * grad_eI) / local_scaling[i][1];
          //double lap_I3 = (dphi[i][qp] * grad_eI3) / local_scaling[i][2];
          //double lap_C = (dphi[i][qp] * grad_eC) / local_scaling[i][3];
          double lap_e = (dphi[i][qp] * grad_en);
          double lap_I = (dphi[i][qp] * grad_eI);
          double lap_I3 = (dphi[i][qp] * grad_eI3);
          double lap_C = (dphi[i][qp] * grad_eC);
          double dsigma_n_x_lap = dsigma_n * lap_e;
          double dsigma_I_x_lap = dsigma_I * lap_I;
          double dsigma_I3_x_lap = dsigma_I3 * lap_I3;
          double dsigma_C_x_lap = dsigma_C * lap_C;

          double volume = elem->volume();

          cons_cat->add(dof_indices_eC[i], -phi0 * J * dC_dphi * phi[i][qp] / C0_C / scaling_C );
          cons_cat->add(dof_indices_u[i], phi0 * J * dC_dphi * phi[i][qp]  / C0_C / scaling_C );

          cons_iodine->add(dof_indices_eI[i], -phi0 * J * dI_dphi * phi[i][qp] / (3.0*C0_tot*scaling_tot) );
          cons_iodine->add(dof_indices_eI3[i], -phi0 * J * dI3_dphi * phi[i][qp] / (C0_tot * scaling_tot) );
          cons_iodine->add(dof_indices_u[i], phi0 * J * phi[i][qp] * ( ( (1/3.0) * dI_dphi + dI3_dphi ) / (C0_tot * scaling_tot) ) );

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
              Knu(i,j) += dsigma_n_x_phi / local_scaling[i][0];
              Knn(i,j) -= dsigma_n_x_phi / local_scaling[i][0];
            }

            if (is_electrolyte && !poisson_only())
            {
              KIu(i,j) += dsigma_I_x_phi / local_scaling[i][1];
              KII(i,j) -= dsigma_I_x_phi / local_scaling[i][1];

              KI3u(i,j) += dsigma_I3_x_phi / local_scaling[i][2];
              KI3I3(i,j) -= dsigma_I3_x_phi / local_scaling[i][2];

              KCu(i,j) += dsigma_C_x_phi / local_scaling[i][3];
              KCC(i,j) -= dsigma_C_x_phi / local_scaling[i][3];
            }




            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            Kuu(i,j) -= drho * phi_i_x_phi_j / local_scaling[i][4];

            if (is_TiO2)
            {
              Kun(i,j) -= phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;

              if (!poisson_only())
              {
                Knu(i,j) -= dR[0] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                Knn(i,j) -= dR[1] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
              }
            }

            if (is_electrolyte)
            {
              KuI(i,j) -= phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuI3(i,j) -= phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuC(i,j) += phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              
              if (!poisson_only())
              {
                if (is_TiO2)
                {
                  KnI(i,j) -= dR[2] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                  KnI3(i,j) -= dR[3] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;

                  KIn(i,j) += 1.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                  KI3n(i,j) -= 0.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                }

                KIu(i,j) += 1.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII(i,j) += 1.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII3(i,j) += 1.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;

                KI3u(i,j) -= 0.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3I(i,j) -= 0.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3I3(i,j) -= 0.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
              }
            }
          }

        }

      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        double J_x_rho = J * sc->get_charge_density() / C0_tot;
        //double J_x_rho = 0.0;

        // net recombination rate
        double J_x_R = 0.0;
        if ((is_TiO2 && is_electrolyte) && !poisson_only())
          J_x_R = J * R;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double net_recomb = J_x_R * phi[i][qp];

          Fu(i) -= J_x_rho * phi[i][qp] / local_scaling[i][4];

          Fn(i) -= net_recomb / R0_e / local_scaling[i][0];

          // TODO factors
          FI(i) += 1.5 * net_recomb / R0_I / local_scaling[i][1];
          FI3(i) -= 0.5 * net_recomb / R0_I3 / local_scaling[i][2];

          double volume = elem->volume();

          tot_cat += J * phi[i][qp] * n_C / C0_C;
          tot_iodine += J * phi[i][qp] * (n_I3 + n_I / 3.0) / C0_tot;

        }
      }

    } // end loop over quadrature points



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
        // we need to know if it is an outer or inner boundary
        bool true_boundary = environment.is_outer_boundary(side);
        bool inner_boundary = environment.is_inner_boundary(side);

        Boundary* boundary = environment.get_boundary(side);

        DSSCContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<DSSCContact*>(
              boundary->get_boundary_properties(get_id()));

        AutoPtr<Elem> side(elem->build_side(s));

        fe_face->reinit(elem, s);



        // calculate the fluxes on the nodes
        if ((contact != NULL) && (dim > 1))
        {
          vector<Point> p(side->n_nodes());
          for (unsigned int i = 0; i < side->n_nodes(); i++)
          {
            nodes_on_boundary_sides.insert(side->node(i));
            p[i] = side->point(i);
          }
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

         if (true_boundary)
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
                u  += phi_face[i][qp] * Xu(i);
                en += phi_face[i][qp] * Xn(i);
                eI += phi_face[i][qp] * XI(i);
                eI3 += phi_face[i][qp] * XI3(i);
                eC += phi_face[i][qp] * XC(i);
                e_field += dphi_face[i][qp] * Xu(i);
                grad_en += dphi_face[i][qp] * Xn(i);
                grad_eI += dphi_face[i][qp] * XI(i);
                grad_eI3 += dphi_face[i][qp] * XI3(i);
                grad_eC += dphi_face[i][qp] * XC(i);
              }

              sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);
              sc->set_coordinates(q_point_face[qp]);

              sc->calculate_densities();
              sc->calculate_traps();
              sc->calculate_equilibrium_traps();

              // we put the phi0 here for convenience
              double sigma_I = sc->get_mobility_I() * sc->get_density_I() / C0_I;
              double dsigma_I = phi0 / C0_I * sc->get_mobility_I() *
                sc->get_density_derivative_I();
              double sigma_I3 = sc->get_mobility_I3() * sc->get_density_I3() / C0_I3;
              double dsigma_I3 = phi0 / C0_I3 * sc->get_mobility_I3() *
                sc->get_density_derivative_I3();
              double sigma_n = sc->get_mobility_n() * sc->get_density_n() / C0_e;
              double dsigma_n = phi0 / C0_e * sc->get_mobility_n() *
                sc->get_density_derivative_n();
              double sigma_C = sc->get_mobility_C() * sc->get_density_C() / C0_C;
              double dsigma_C = phi0 / C0_C * sc->get_mobility_C() *
                sc->get_density_derivative_C();


            //define flux direction
            

            // the jacobian x weight x scaling
            double J = JxW_face[qp];

            if (contact->is_cathode())
            {
              for (unsigned int i = 0; i < n_dofs; i++)
              {

                if (jacobian != NULL)
                {
                  for (unsigned int j = 0; j < n_dofs; j++)
                  {

                    //double res = contact->get_load() * x0;
                    double j0 = contact->get_ex_curr();
                
                    double kT = sc->get_lattice_temperature();
                    double kT2 = 2*kT;
	            //double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	            //double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
                    double Normal_I = x0 / (phi0 * C0_I * Constants::e * scaling_I);
                    double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * scaling_I3);
                
                    //double Eredox = -2*u + 1.5*eI - 0.5*eI3;
                    double Eredox = -u;

                    double J_phi_i_phi_j = J * phi_face[i][qp] * phi_face[j][qp];
	        
                    //KIu(i,j) += 1.5 * 2 * J_phi_i_phi_j * Normal_I / res;
                    //KII(i,j) += 1.5 * -1.5 * J_phi_i_phi_j * Normal_I / res;
                    //KII3(i,j) += 1.5 * 0.5 * J_phi_i_phi_j * Normal_I / res;
//                    KIu(i,j) += 1.5 * J_phi_i_phi_j * Normal_I / res;
                    KIu(i,j) += 1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I / kT2;

                    //KI3u(i,j) += -0.5 * 2 * J_phi_i_phi_j * Normal_I3 / res;
                    //KI3I(i,j) += -0.5 * -1.5 * J_phi_i_phi_j * Normal_I3 / res;
                    //KI3I3(i,j) += -0.5 * 0.5 * J_phi_i_phi_j * Normal_I3 / res;
//                    KI3u(i,j) += -0.5 * J_phi_i_phi_j * Normal_I3 / res;
                    KI3u(i,j) += -0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I3 / kT2;

                  }
                }
                if (residual != NULL)
                {
                    //double res = contact->get_load() * x0;
                    double j0 = contact->get_ex_curr();
	            //double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	            //double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
                    double Normal_I = x0 / (phi0 * C0_I * Constants::e * scaling_I);
                    double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * scaling_I3);
                    double kT = sc->get_lattice_temperature();
                    double kT2 = 2*kT;
               
                    //double Eredox = -2*u + 1.5*eI - 0.5*eI3;
                    double Eredox = -u;

                    FI(i) += -1.5 * j0 * ( exp(Eredox/kT2)   -  exp(-Eredox/kT2) ) * J * phi_face[i][qp] * Normal_I;
                    FI3(i) += 0.5 * j0 * ( exp(Eredox/kT2)   -  exp(-Eredox/kT2) ) * J * phi_face[i][qp] * Normal_I3;

                  //pot = -u;
                  //pot = -2*u + 1.5*eI - 0.5*eI3;
//                  FI(i) += -1.5 * pot * J * phi_face[i][qp] * Normal_I / res ;
//                  FI3(i) += 0.5 * pot * J * phi_face[i][qp] * Normal_I3 / res ;
                }
              }

            }
            else
            {

             if (n_cat == NULL)
             {
               n_cat = side->get_node(0);
             }

            for (unsigned int i = 0; i < n_dofs; i++)
            {
               if (jacobian != NULL)
               {
               for (unsigned int j = 0; j < n_dofs; j++)
                {
	           //double Normal_n = x0 / (phi0 * C0_e * local_scaling[s][0] );
	           double Normal_n = x0 / (phi0 * C0_e );
                   double kinetic_anode = contact->get_kinetic();
                    
                   double J_phi_i_phi_j = J * phi_face[i][qp] * phi_face[j][qp];

                   double dn_dphi = sc->get_density_derivative_n();

                   Knu(i,j) += -kinetic_anode * dn_dphi * Normal_n * J_phi_i_phi_j;
                   Knn(i,j) += -kinetic_anode * -dn_dphi * Normal_n * J_phi_i_phi_j;
                }
               }
               if (residual != NULL)
               { 
                 
                 double density_n = sc->get_density_n();
                 double n_dark = sc->get_equilibrium_concentrations().n;
                 double bias = contact->get_potential();
                 double kT = sc->get_lattice_temperature();
                 double n_0 = n_dark * exp(bias/kT); 

	         //double Normal_n = x0 / ( phi0 * C0_e * local_scaling[s][0] );
	         double Normal_n = x0 / ( phi0 * C0_e );
                 double kinetic_anode = contact->get_kinetic();

                 Fn(i) += -kinetic_anode * (density_n - n_0) * Normal_n * J * phi_face[i][qp];
               }
             }

            }

         } //contact

        } // qp

        } // true_boundary
        } // dimension
        else // i.e. dim == 1
        {


        if (true_boundary)
        {
          nodes_on_boundary_sides.insert(elem->node(s));

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
            grad_en += dphi_face[n][0](0) * Xn(n);
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
          sc->calculate_traps();
          sc->calculate_equilibrium_traps();

          // we put the phi0 here for convenience
          double sigma_I = sc->get_mobility_I() * sc->get_density_I() / C0_I;
          double dsigma_I = phi0 / C0_I * sc->get_mobility_I() *
            sc->get_density_derivative_I();
          double sigma_I3 = sc->get_mobility_I3() * sc->get_density_I3() / C0_I3;
          double dsigma_I3 = phi0 / C0_I3 * sc->get_mobility_I3() *
            sc->get_density_derivative_I3();
          double sigma_n = sc->get_mobility_n() * sc->get_density_n() / C0_e;
          double dsigma_n = phi0 / C0_e * sc->get_mobility_n() *
            sc->get_density_derivative_n();
          double sigma_C = sc->get_mobility_C() * sc->get_density_C() / C0_C;
          double dsigma_C = phi0 / C0_C * sc->get_mobility_C() *
            sc->get_density_derivative_C();

          double x_c = elem->centroid()(0);
          double x_s = elem->point(s)(0);
          double sign = (x_s > x_c) ? 1 : -1;


          if (contact->is_cathode())
          {
//             if (n_cat == NULL)
//             {
//               n_cat = side->get_node(0);            
//             }

            if (jacobian != NULL)
             {

                //double res = contact->get_load() * x0;
                
                double j0 = contact->get_ex_curr();
                double kT = sc->get_lattice_temperature();
                double kT2 = 2*kT;
	        double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	        double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
                
                double Eredox = -Xu(s);
               
                double bias = contact->get_potential();
                
                double Eredox1 = bias + 2.5*XI(s) - 1.5*XI3(s) - 2*Xu(s);
                double Eredox2 = -bias +4*Xu(s) - 3.5*XI(s) + 0.5*XI3(s);

//                KI3u(s,s) += 0.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//                KI3I(s,s) += 0.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//                KI3I3(s,s) += 0.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
                
//                KIu(s,s) += -1.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//                KII(s,s) += -1.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//                KII3(s,s) += -1.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
	        

//               KIu(s,s) += 1.5 * Normal_I / res;
//               KI3u(s,s) += -0.5 * Normal_I3 / res;
                
                KIu(s,s) += -1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I / kT2;

                KI3u(s,s) += 0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I3 / kT2;
                

            }
            if (residual != NULL)
            {
//               double res = contact->get_load() * x0;

	       double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	       double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
               double kT = sc->get_lattice_temperature();
               double kT2 = 2*kT;
               
               double j0 = contact->get_ex_curr();
               double Eredox = -Xu(s);
               
               double bias = contact->get_potential();
               
//               double Eredox1 = bias + 2.5*XI(s) - 1.5*XI3(s) - 2*Xu(s);
//               double Eredox2 = -bias +4*Xu(s) - 3.5*XI(s) + 0.5*XI3(s);
               //pot =  -Xu(s);
//               FI(s) += -1.5 * pot * Normal_I / res ;
//               FI3(s) += 0.5 * pot * Normal_I3 / res ;

               FI3(s) += 0.5 * j0 * ( exp(Eredox/kT2)   -  exp(-Eredox/kT2) ) * Normal_I3;
               FI(s) += -1.5 * j0 * ( exp(Eredox/kT2)   -  exp(-Eredox/kT2) ) * Normal_I;
//               FI3(s) += 0.5 * j0 * ( exp(Eredox1/kT2)   -  exp(Eredox2/kT2) ) * Normal_I3;
//               FI(s) += -1.5 * j0 * ( exp(Eredox1/kT2)   -  exp(Eredox2/kT2) ) * Normal_I;

             }
           }
           else
	   {

             if (n_cat == NULL)
             {
               n_cat = side->get_node(0);            
             }

             if (jacobian != NULL)
             {
               
	       //double Normal_n = x0 / (phi0 * C0_e * Constants::e * local_scaling[s][0] );
	       double Normal_n = x0 / (phi0 * C0_e * local_scaling[s][0] );
               double kinetic_anode = contact->get_kinetic();
               
               double bias = contact->get_potential();
               double kT = sc->get_lattice_temperature();
               double j0 = contact->get_ex_curr();
               double kT2 = 2*kT;
	       double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	       double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
               
               double Eredox1 = bias + 2.5*XI(s) - 1.5*XI3(s) - 2*Xu(s);
               double Eredox2 = -bias +4*Xu(s) - 3.5*XI(s) + 0.5*XI3(s);
                 
               double dn_dphi = sc->get_density_derivative_n();
               
//               KI3u(s,s) += 0.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//               KI3I(s,s) += 0.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//               KI3I3(s,s) += 0.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
               
//               KIu(s,s) += -1.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//               KII(s,s) += -1.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//               KII3(s,s) += -1.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;

               Knu(s,s) += -kinetic_anode * dn_dphi * Normal_n;
               Knn(s,s) += -kinetic_anode * -dn_dphi * Normal_n;

             }
             if (residual != NULL)
             { 
                 
               double density_n = sc->get_density_n();
               double n_dark = sc->get_equilibrium_concentrations().n;
               double bias = contact->get_potential();
               double kT = sc->get_lattice_temperature();
               double n_0 = n_dark * exp(bias/kT); 
               
               double j0 = contact->get_ex_curr();
               double kT2 = 2*kT;
	       double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	       double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
               
               double Eredox1 = bias + 2.5*XI(s) - 1.5*XI3(s) - 2*Xu(s);
               double Eredox2 = -bias +4*Xu(s) - 3.5*XI(s) + 0.5*XI3(s);

	       //double Normal_n = x0 / ( phi0 * C0_e * Constants::e * local_scaling[s][0] );
	       double Normal_n = x0 / ( phi0 * C0_e * local_scaling[s][0] );
               double kinetic_anode = contact->get_kinetic();

               Fn(s) += -kinetic_anode * (density_n - n_0) * Normal_n; 
               
//               FI3(s) += 0.5 * j0 * ( exp(Eredox1/kT2)   -  exp(Eredox2/kT2) ) * Normal_I3;
//               FI(s) += -1.5 * j0 * ( exp(Eredox1/kT2)   -  exp(Eredox2/kT2) ) * Normal_I;

             }

           }
	}
       } // if (dim ...)
      }
    } // end loop over element sides


    // check if it is a dielectric
    if (!sc->is_TiO2())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          Kun(i, j) = 0.0;
          Knu(i, j) = Knn(i, j) = KnI(i, j) = KnI3(i, j) = 0.0;
          KIn(i, j) = KI3n(i, j) = 0.0;
        }

        if (!is_internal_boundary_node(elem->get_node(i)))
          Knn(i, i) = 1.0;

        // we simply set the electrochemical potential to zero
        Fn(i) = 0.0;
      }
    }

    // check if it is a dielectric
    if (!sc->is_electrolyte())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          KuI(i, j) = KuI3(i, j) = KuC(i, j) = 0.0;

          KIu(i, j) = KIn(i, j) = KII(i, j) = KII3(i, j) = 0.0;
          KI3u(i, j) = KI3n(i, j) = KI3I(i, j) = KI3I3(i, j) = 0.0;
          KCu(i, j) = KCC(i, j) = 0.0;
          
          KnI(i, j) = KnI3(i, j) = 0.0;
        }

        if (!is_internal_boundary_node(elem->get_node(i)))
        {
          KII(i, i) = 1.0;
          KI3I3(i, i) = 1.0;
          KCC(i, i) = 1.0;
        }

        // we simply set the electrochemical potential to zero
        FI(i) = 0.0;
        FI3(i) = 0.0;
        FC(i) = 0.0;
      }
    } 

    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    // NOTE: this changes dof_indices that's why th application of
    //       Dirichlet type BCs needs special care
    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

///*
    //
    // now as last thing we apply Dirichlet type Bcs
    //
    BoundaryNodeList::const_iterator node_it;
    const BoundaryNodeList::const_iterator end =
      dirichlet_nodes.end();
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

          // we only impose Dirichlet type BCs if the node has an associated
          // boundary side
          if (nodes_on_boundary_sides.find(elem->node(i)) !=
              nodes_on_boundary_sides.end())
          {
            if ( (!contact->is_cathode()) && (!contact->is_gate()) )
            {

              double bias = contact->get_potential();

//              Ke.condense(i + n_dofs, i + n_dofs, bias, Fe);

            }
            else if (contact->is_gate())
            {
               double bias = contact->get_potential();

              Ke.condense(i, i, bias, Fe);
            }
          }
        }
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


  unsigned int dof_cat = n_cat->dof_number(system.number(), eC_var, 0);
  unsigned int dof_iodine = n_cat->dof_number(system.number(), eI_var, 0);

  if (jacobian != NULL)
  {

    jacobian->close();

    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_cat->size(); j++)
      jacobian->set(dof_cat, j, (*cons_cat)(j));

    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_iodine->size(); j++)
      jacobian->set(dof_iodine, j, (*cons_iodine)(j));

    jacobian->close();
//    jacobian->print_matlab("J.m");
  }
  else
  {
    residual->close();
    residual->set(dof_cat, (tot_cat / scaling_C - _cation_amount / C0_C / scaling_C));
    residual->set(dof_iodine, (tot_iodine / scaling_tot - _iodine_amount / C0_tot / scaling_tot));
    residual->close();
//    residual->print_matlab("F.m");
  }


  perf_log.stop_event("assembly");

}





void
DSSC::do_setup_solution_variables(void)
{

  if (plot_solution("Potential"))
  { 
     add_plot_variable("ElPotential");
     add_plot_variable("eQFermi");
     add_plot_variable("IQFermi");
     add_plot_variable("I3QFermi");
     add_plot_variable("CQFermi");
     add_plot_variable("Eredox");
  }

  if (plot_solution("Density"))
  { 
     add_plot_variable("eDensity");
     add_plot_variable("IDensity");
     add_plot_variable("I3Density");
     add_plot_variable("CDensity");
     add_plot_variable("Traps");
     add_plot_variable("NetRecombination");
     add_plot_variable("Generation");
  }

  if (plot_solution("Current"))
  { 
     add_plot_variable("ElField");
     add_plot_variable("CurrentDensity");
     add_plot_variable("eCurrentDensity");
     add_plot_variable("ICurrentDensity");
     add_plot_variable("I3CurrentDensity");
     add_plot_variable("CCurrentDensity");
  }  

  if (plot_solution("Mobility"))
  {
     add_plot_variable("eMobility");
  }


  // the contact currents/voltages
   
  // declare solution variables
  declare_solution(ElPotential, REAL, NODES, "V");
  declare_solution(eQFermi, REAL, NODES, "eV");
  declare_solution(IQFermi, REAL, NODES, "eV");
  declare_solution(I3QFermi, REAL, NODES, "eV");
  declare_solution(CQFermi, REAL, NODES, "eV");
  declare_solution(Eredox, REAL, NODES, "eV");
  
  declare_solution(ElField, VECTOR, CELL, "V/cm");

  declare_solution(eDensity, REAL, NODES, "cm^-3");
  declare_solution(IDensity, REAL, NODES, "cm^-3");
  declare_solution(I3Density, REAL, NODES, "cm^-3");
  declare_solution(CDensity, REAL, NODES, "cm^-3");
  declare_solution(Traps, REAL, NODES, "cm^-3");

  //declare_solution(eMobility, REAL, NODES, "cm^2/(V*s)");

  //declare_solution(eConductivity, REAL, NODES, "S/cm");

  declare_solution(CurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(eCurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(ICurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(I3CurrentDensity, VECTOR, CELL, "A/cm^2");
  declare_solution(CCurrentDensity, VECTOR, CELL, "A/cm^2");

  //declare_solution(IonizedElectronTraps, REAL, NODES, "cm^-3");

  declare_solution(Generation, REAL, NODES, "1/(s*cm^3)");
  declare_solution(NetRecombination, REAL, NODES, "1/(s*cm^3)");
  
  declare_solution(eMobility, REAL, NODES, "cm^2/(V*s)");

    //bool plot_rec = plot_solution(NetRecombination);
/*
    const set<PhysicalModel*>& pm = get_physical_models();
    set<PhysicalModel*>::const_iterator it(pm.begin());
    set<PhysicalModel*>::const_iterator end(pm.end());

    for ( ; it != end; ++it)
    {
      DSSCModel* sc =
          static_cast<DSSCModel*>(*it);

     


      //vector<ID> ids;
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
*/
  
  // the contact currents/voltages
   
/*    unsigned int dim = get_mesh().mesh_dimension();
    
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
*/
//    declare_solution(ContactCurrent, REAL, GLOBAL, units);
//    add_alias("ContactCurrents", ContactCurrent);
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

        //const DSSCContact* mod = static_cast<const DSSCContact*>(modit->second);
        //const DSSCContact* mod =
        //  static_cast<const DSSCContact*>(it->first->get_boundary_properties(get_id()));
        //  const DSSCContact* mod =
         //  static_cast<const DSSCContact*>(modit->second);
        //if (mod->has_current())
        //{
          string name(bd->get_name() + ".current");
          // if currents should be plotted, add it also to the plot variables
          if (plot_curr) add_plot_variable(name);
          cnt_set.insert(name);
        //}
      }
    }
    string units = "A/cm^2";
   
    // now we declare them
    unsigned int id = static_cast<ID>(ContactCurrent);
    for (set<string>::iterator i(cnt_set.begin()); i != cnt_set.end(); ++i)
    {
      ++id;
      declare_solution_ext(*i, id, SolutionDescriptor::REAL,
          SolutionDescriptor::GLOBAL, units);
    }

}







void
DSSC::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& points)
{

  unsigned int np = points.size();

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const NumericVector<Number>& solution = system->get_solution_vector();
//  const NumericVector<Number>& oldsolution = system->get_vector("old_sol");

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int u_var = system->variable_number("potential");
  unsigned int en_var = system->variable_number("fermi_n");
  unsigned int eI_var = system->variable_number("fermi_I");
  unsigned int eI3_var = system->variable_number("fermi_I3");
  unsigned int eC_var = system->variable_number("fermi_C");

  FEType fe_type = system->variable_type(u_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_u;
  vector<unsigned int> dof_indices_en;
  vector<unsigned int> dof_indices_eI;
  vector<unsigned int> dof_indices_eI3;
  vector<unsigned int> dof_indices_eC;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  DSSCModel* sc =
    dynamic_cast<DSSCModel*>(
        get_physical_model(subdomain));

  assert(sc != NULL);

  sc->reinit(elem);

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_u, u_var);
  dof_map.dof_indices(elem, dof_indices_en, en_var);
  dof_map.dof_indices(elem, dof_indices_eI, eI_var);
  dof_map.dof_indices(elem, dof_indices_eI3, eI3_var);
  dof_map.dof_indices(elem, dof_indices_eC, eC_var);

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
  RealGradient jI(0);
  RealGradient jI3(0);
  RealGradient jC(0);
  RealGradient el_field(0);

  for (unsigned int n = 0; n < np; n++)
  {
    double u  = 0.0;
    double en = 0.0;
    double eI = 0.0;
    double eI3 = 0.0;
    double eC = 0.0;
    RealGradient e_field(0);
    RealGradient grad_en_loc(0);
    RealGradient grad_eI_loc(0);
    RealGradient grad_eI3_loc(0);
    RealGradient grad_eC_loc(0);

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      u += phi[i][n] * solution(dof_indices_u[i]);
      en += phi[i][n] * solution(dof_indices_en[i]);
      eI += phi[i][n] * solution(dof_indices_eI[i]);
      eI3 += phi[i][n] * solution(dof_indices_eI3[i]);
      eC += phi[i][n] * solution(dof_indices_eC[i]);

      grad_en_loc += dphi[i][n] * solution(dof_indices_en[i]);
      grad_eI_loc += dphi[i][n] * solution(dof_indices_eI[i]);
      grad_eI3_loc += dphi[i][n] * solution(dof_indices_eI3[i]);
      grad_eC_loc += dphi[i][n] * solution(dof_indices_eC[i]);

      e_field += dphi[i][n] * solution(dof_indices_u[i]);

    }

    // scale the potential back
    u *= phi0;
    en *= phi0;
    eI *= phi0;
    eI3 *= phi0;
    eC *= phi0;
    e_field *= -phi0;
    grad_en_loc *= phi0;
    grad_eI_loc *= phi0;
    grad_eI3_loc *= phi0;
    grad_eC_loc *= phi0;

    //double kT = sc->get_lattice_temperature();

    sc->set_coordinates(real_pts[n]);

    sc->set_potentials(u, en, eI, eI3, eC);

    sc->set_electric_field(e_field);
    sc->set_grad_fermi_n(grad_en_loc);
    sc->set_grad_fermi_I(grad_eI_loc);
    sc->set_grad_fermi_I3(grad_eI3_loc);
    sc->set_grad_fermi_C(grad_eC_loc);

    sc->calculate_densities();
    sc->calculate_traps();
    sc->calculate_equilibrium_traps();
    sc->calculate_net_recombination_rate();

    double edens = sc->get_density_n_exp_DOS();
    double edens_c = sc->get_density_n();
    double Idens = sc->get_density_I();
    double I3dens = sc->get_density_I3();
    double Cdens = sc->get_density_C();

    sc->calculate_net_recombination_rate();

    double sigma_e = Constants::e * edens_c * sc->get_mobility_n();
    double sigma_I = Constants::e * Idens * sc->get_mobility_I();
    double sigma_I3 = Constants::e * I3dens * sc->get_mobility_I3();
    double sigma_C = Constants::e * Cdens * sc->get_mobility_C();


    RealGradient dfn = grad_en_loc;
    RealGradient dfI = grad_eI_loc;
    RealGradient dfI3 = grad_eI3_loc;
    RealGradient dfC = grad_eC_loc;
    RealGradient jn_loc = -sigma_e * dfn;
    RealGradient jI_loc = -sigma_I * dfI;
    RealGradient jI3_loc = -sigma_I3 * dfI3;
    RealGradient jC_loc = sigma_C * dfC;
    jn += jn_loc;
    jI += jI_loc;
    jI3 += jI3_loc;
    jC += jC_loc;

    el_field += e_field;


    if (values.count(ElPotential))
      values[ElPotential][n] = u;

    if (values.count(eQFermi))
      values[eQFermi][n] = -en;

    if (values.count(IQFermi))
      values[IQFermi][n] = -eI;

    if (values.count(I3QFermi))
      values[I3QFermi][n] = -eI3;

    if (values.count(CQFermi))
      values[CQFermi][n] = -eC;

    if (values.count(Eredox))
      values[Eredox][n] = 0.5*(3*eI - eI3 - 2*u);

//    if (values.count(Ec))
//      values[Ec][n] = sc->get_conduction_band_edge() - u;

//    if (values.count(Ev))
//      values[Ev][n] = sc->get_valence_band_edge() - u;

//    if (values.count(Ec0))
//      values[Ec0][n] = sc->get_conduction_band_edge();

//    if (values.count(Ev0))
//      values[Ev0][n] = sc->get_valence_band_edge();

//    if (values.count(Eg))
//      values[Eg][n] =
//        sc->get_conduction_band_edge() - sc->get_valence_band_edge();

    if (values.count(eDensity))
      values[eDensity][n] = edens;

    if (values.count(IDensity))
      values[IDensity][n] = Idens;

    if (values.count(I3Density))
      values[I3Density][n] = I3dens;
    
    if (values.count(CDensity))
      values[CDensity][n] = Cdens;
    
    if (values.count(Traps))
      values[Traps][n] = - sc->get_ionized_electron_traps();
    
    //if (values.count(eMobility))
    //  values[eMobility][n] = sc->get_electron_mobility();


    //if (values.count(eConductivity))
    //  values[eConductivity][n] = sigma_e;

    //if (values.count(hConductivity))
    //  values[hConductivity][n] = sigma_h;


    //bool trapped_electrons = values.count(IonizedElectronTraps);
    //if (trapped_electrons || trapped_holes)
    //{
    //  sc->calculate_traps();
    //  if (trapped_electrons)
    //    values[IonizedElectronTraps][n] = sc->get_ionized_electron_traps();
    //}

    
     // bool get_recomb = values.count(NetRecombination);
     // double tot_rec = 0;

      // loop over all recombination ids
      //set<ID>::const_iterator rec_it(_recombination_ids.begin());
      //for ( ; rec_it != _recombination_ids.end(); ++rec_it)
      //{
      //  bool requested = values.count(*rec_it);
      //  double rec = 0;
      //  if (get_recomb || requested)
      //  rec = sc->get_net_recombination_rate();
      // 
      //  if (requested)
      //    values[*rec_it][n] = rec;

    //if (get_recomb)
    //  tot_rec += sc->get_recombination_rate();  

    if (values.count(NetRecombination))
      values[NetRecombination][n] = sc->get_recombination_rate();
    
    if (values.count(Generation))
      values[Generation][n] = sc->get_generation_rate();

    if (values.count(eMobility))
      values[eMobility][n] = sc->get_dens_elec_mobility();
  }


  // the cell based solutions
  // for now, we make a mean value

 
  if (values.count(ElField))
  {
    values[ElField][0] = el_field(0) / np;
    values[ElField][1] = el_field(1) / np;
    values[ElField][2] = el_field(2) / np;
  }

  if (values.count(CurrentDensity))
  {
    values[CurrentDensity][0] = (jn(0) + jI(0) + jI3(0)) / np;
    values[CurrentDensity][1] = (jn(1) + jI(1) + jI3(1)) / np;
    values[CurrentDensity][2] = (jn(2) + jI(2) + jI3(2)) / np;
  }

  if (values.count(eCurrentDensity))
  {
    values[eCurrentDensity][0] = jn(0) / np;
    values[eCurrentDensity][1] = jn(1) / np;
    values[eCurrentDensity][2] = jn(2) / np;
  }

  if (values.count(ICurrentDensity))
  {
    values[ICurrentDensity][0] = jI(0) / np;
    values[ICurrentDensity][1] = jI(1) / np;
    values[ICurrentDensity][2] = jI(2) / np;
  }

  if (values.count(I3CurrentDensity))
  {
    values[I3CurrentDensity][0] = jI3(0) / np;
    values[I3CurrentDensity][1] = jI3(1) / np;
    values[I3CurrentDensity][2] = jI3(2) / np;
  }

  if (values.count(CCurrentDensity))
  {
    values[CCurrentDensity][0] = jC(0) / np;
    values[CCurrentDensity][1] = jC(1) / np;
    values[CCurrentDensity][2] = jC(2) / np;
  }
  
//  if (values.count(ConductionBands))
//  {
//    const vector<double>& cb = sc->get_conduction_bands();
//    for (size_t i = 0; i < cb.size(); ++i)
//      values[ConductionBands][i] = cb[i];
//  }

}




void
DSSC::get_solution_secure(map<ID, vector<double> >& values)
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
DSSC::do_assembly_frequency(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  PerfLog perf_log("Matrix assembly", false);
  perf_log.start_event("assembly EIS");

  // references for nicer code
  const MeshBase& mesh = get_mesh();
  TiberNonlinearSystem& system_frequency = get_equation_system<TiberNonlinearSystem>(1);

  /// COMMENT: here we define a dumb value for freq (frequency for the EIS)
  double freq = 0.0;

  // references for nicer code
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);
  //get_solution (riprendere la parte della soluzione)
  
  // finally initialize the newly created system
  system_frequency.init();

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
  double C0_e = _cond_scaling.n;
  //double C0_e = _cond_scaling.C;
  // scaling for I
  double C0_I = _cond_scaling.I;
  // scaling for I3
  double C0_I3 = _cond_scaling.I3;
  // scaling for C
  double C0_C = _cond_scaling.C;
  // scaling for I + I3
  double C0_tot = C0;
  // scaling C
  double scaling_C = 1.0;
  // scaling I
  double scaling_I = 1.0;
  // scaling I3
  double scaling_I3 = 1.0;
  // scaling all charges
  double scaling_tot = 1.0;



  if (do_local_scaling_)
  {
    scaling_C = _cond_scaling.C;
    scaling_I = _cond_scaling.I;
    scaling_I3 = _cond_scaling.I3;
    scaling_tot = scaling.get_density_scaling();
    C0_tot = C0_e = C0_I = C0_I3 = C0_C = 1.0;
  }



  // scaling for recombination rates
  double tmp = phi0 / (x0 * x0);
  double R0_e = C0_e * tmp;
  double R0_I = C0_I * tmp ;
  double R0_I3 = C0_I3 * tmp;


  const DofMap& dof_map = system_frequency.get_dof_map();

  // numeric ids corresponding to the variables
  const unsigned int u_var_R = system.variable_number("potential_R");
  const unsigned int en_var_R = system.variable_number("fermi_n_R");
  const unsigned int eI_var_R = system.variable_number("fermi_I_R");
  const unsigned int eI3_var_R = system.variable_number("fermi_I3_R");
  const unsigned int eC_var_R = system.variable_number("fermi_C_R");
  const unsigned int u_var_I = system.variable_number("potential_I");
  const unsigned int en_var_I = system.variable_number("fermi_n_I");
  const unsigned int eI_var_I = system.variable_number("fermi_I_I");
  const unsigned int eI3_var_I = system.variable_number("fermi_I3_I");
  const unsigned int eC_var_I = system.variable_number("fermi_C_I");

  FEType fe_type = system_frequency.variable_type(u_var_R);

  libMeshEnums::Order integration_order = libMeshEnums::FIFTH;
  //libMeshEnums::Order integration_order = libMeshEnums::FIRST;

  // the finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, integration_order);
  //QTrap qrule(dim);
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
  // the steady-state solution
  DenseVector<Number> S;
    
    // Reposition the submatrices according to this scheme:
    //
    //                                                                                     
    //        | KuuR  KunR  KuIR  KuI3R  KuCR   0      0      0      0        0     |        | 0   |
    //        | KnuR  KnnR  KnIR  KnI3R   0    KnuFR  KnnFR   0      0        0     |        | FnR |
    //        | KIuR  KInR  KIIR  KII3R   0    KIuFR   0     KIIFR   0        0     |        | 0   |
    //        | KI3uR KI3nR KI3IR KI3I3R  0    KI3uFR  0      0     KI3I3FR   0     |        | 0   |
    //   Ke = | KCuR   0     0     0     KCCR  KCuFR   0      0      0       KCCFR  |;  Fe = | 0   |
    //        |  0     0     0     0      0    KuuI   KunI   KuII   KuI3I    KuCI   |        | 0   |
    //        | KnuFI KnnFI  0     0      0    KnuI   KnnI   KnII   KnI3I     0     |        | FnI |
    //        | KIuFI  0    KIIFI  0      0    KIuI   KInI   KIII   KII3I     0     |        | 0   |
    //        | KI3uFI 0     0    KI3I3FI 0    KI3uI  KI3nI  KI3II  KI3I3I    0     |        | 0   |
    //        | KCuFI   0    0     0     KCCFI KCuI    0      0      0       KCCI   |        | 0   |
    //                               
    //

  DenseSubMatrix<Number>
    KuuR(Ke), KunR(Ke), KuIR(Ke), KuI3R(Ke), KuCR(Ke), /* KFuuR = 0, KFunR = 0, KFuIR = 0, KFuI3R = 0, KFuCR = 0, */ 
    KnuR(Ke), KnnR(Ke), KnIR(Ke), KnI3R(Ke), /* KnCR = 0, */ KnuFR(Ke), KnnFR(Ke), /* =0, =0, =0, */
    KIuR(Ke), KInR(Ke), KIIR(Ke), KII3R(Ke), /* KICR = 0, */ KIuFR(Ke),  /* =0, */  KIIFR(Ke), /* =0, =0, */
    KI3uR(Ke),KI3nR(Ke),KI3IR(Ke),KI3I3R(Ke), /* KI3CR = 0, */ KI3uFR(Ke), /* =0, =0, */ KI3I3FR(Ke), /* =0, */
    KCuR(Ke), /* =0, =0, =0, */  KCCR(Ke), KCuFR(Ke), /* =0, =0, =0, */ KCCFR(Ke),
    /* =0, =0, =0, =0, =0, */ KuuI(Ke), KunI(Ke), KuII(Ke), KuI3I(Ke), KuCI(Ke),  
    KnuFI(Ke), KnnFI(Ke), /* =0, =0, =0, */ KnuI(Ke), KnnI(Ke), KnII(Ke), KnI3I(Ke), /* =0, */ 
    KIuFI(Ke), /* =0, */ KIIFI(Ke), /* =0, =0, */ KIuI(Ke), KInI(Ke), KIII(Ke), KII3I(Ke), /* =0, */ 
    KI3uFI(Ke), /* =0, =0, */ KI3I3FI(Ke), /* =0, */ KI3uI(Ke), KI3nI(Ke), KI3II(Ke), KI3I3I(Ke), /* =0, */ 
    KCuFI(Ke), /* =0, =0, =0, */ KCCFI(Ke), KCuI(Ke), /* =0,  =0,  =0, */  KCCI(Ke);

  DenseSubVector<Number>
    FuR(Fe),
    FnR(Fe),
    FIR(Fe),
    FI3R(Fe),
    FCR(Fe),
    FuI(Fe),
    FnI(Fe),
    FII(Fe),
    FI3I(Fe),
    FCI(Fe);

  DenseSubVector<Number>
    XuR(X),
    XnR(X),
    XIR(X),
    XI3R(X),
    XCR(X),
    XuI(X),
    XnI(X),
    XII(X),
    XI3I(X),
    XCI(X);
  
  DenseSubVector<Number>
    Su(X),
    Sn(X),
    SI(X),
    SI3(X),
    SC(X);


  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_u_R;
  vector<unsigned int> dof_indices_en_R;
  vector<unsigned int> dof_indices_eI_R;
  vector<unsigned int> dof_indices_eI3_R;
  vector<unsigned int> dof_indices_eC_R;
  vector<unsigned int> dof_indices_u_I;
  vector<unsigned int> dof_indices_en_I;
  vector<unsigned int> dof_indices_eI_I;
  vector<unsigned int> dof_indices_eI3_I;
  vector<unsigned int> dof_indices_eC_I;

  // zero out residual and jacobian !! IMPORTANT !!
  if (residual != NULL)
    residual->zero();
  if (jacobian != NULL)
    jacobian->zero();


  const Node* n_cat = NULL;
  double tot_cat = 0;
  AutoPtr<NumericVector<Number> > cons_cat = x.clone();
  cons_cat->zero();

  const Node* n_iodine = NULL;
  double tot_iodine = 0;
  AutoPtr<NumericVector<Number> > cons_iodine = x.clone();
  cons_iodine->zero();

  set<unsigned int> inner_boundary_nodes;

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
    dof_map.dof_indices(elem, dof_indices_u_R, u_var_R);
    dof_map.dof_indices(elem, dof_indices_en_R, en_var_R);
    dof_map.dof_indices(elem, dof_indices_eI_R, eI_var_R);
    dof_map.dof_indices(elem, dof_indices_eI3_R, eI3_var_R);
    dof_map.dof_indices(elem, dof_indices_eC_R, eC_var_R);
    dof_map.dof_indices(elem, dof_indices_u_I, u_var_I);
    dof_map.dof_indices(elem, dof_indices_en_I, en_var_I);
    dof_map.dof_indices(elem, dof_indices_eI_I, eI_var_I);
    dof_map.dof_indices(elem, dof_indices_eI3_I, eI3_var_I);
    dof_map.dof_indices(elem, dof_indices_eC_I, eC_var_I);

    unsigned int n_dofs     = dof_indices_u_R.size();
    unsigned int n_dofs_tot = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs_tot, n_dofs_tot);
    Fe.resize(n_dofs_tot);
    X.resize(n_dofs_tot);

    // extract local solution, accounting for constraints
    dof_map.extract_local_vector(x, dof_indices, X);

    // Reposition the submatrices according to this scheme:
    //
    //                                                                                     
    //        | KuuR  KunR  KuIR  KuI3R  KuCR   0      0      0      0        0     |        | 0   |
    //        | KnuR  KnnR  KnIR  KnI3R   0    KnuFR  KnnFR   0      0        0     |        | FnR |
    //        | KIuR  KInR  KIIR  KII3R   0    KIuFR   0     KIIFR   0        0     |        | 0   |
    //        | KI3uR KI3nR KI3IR KI3I3R  0    KI3uFR  0      0     KI3I3FR   0     |        | 0   |
    //   Ke = | KCuR   0     0     0     KCCR  KCuFR   0      0      0       KCCFR  |;  Fe = | 0   |
    //        |  0     0     0     0      0    KuuI   KunI   KuII   KuI3I    KuCI   |        | 0   |
    //        | KnuFI KnnFI  0     0      0    KnuI   KnnI   KnII   KnI3I     0     |        | FnI |
    //        | KIuFI  0    KIIFI  0      0    KIuI   KInI   KIII   KII3I     0     |        | 0   |
    //        | KI3uFI 0     0    KI3I3FI 0    KI3uI  KI3nI  KI3II  KI3I3I    0     |        | 0   |
    //        | KCuFI   0    0     0     KCCFI KCuI    0      0      0       KCCI   |        | 0   |
    //                               
    //
    KuuR.reposition(0, 0, n_dofs, n_dofs);
    KunR.reposition(0, n_dofs, n_dofs, n_dofs);
    KuIR.reposition(0, 2 * n_dofs, n_dofs, n_dofs);
    KuI3R.reposition(0, 3 * n_dofs, n_dofs, n_dofs);
    KuCR.reposition(0, 4 * n_dofs, n_dofs, n_dofs);
    //
    KnuR.reposition(n_dofs, 0, n_dofs, n_dofs);
    KnnR.reposition(n_dofs, n_dofs, n_dofs, n_dofs);
    KnIR.reposition(n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KnI3R.reposition(n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    KnuFR.reposition(n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KnnFR.reposition(n_dofs, 6 * n_dofs, n_dofs, n_dofs);
    //
    KIuR.reposition(2 * n_dofs, 0, n_dofs, n_dofs);
    KInR.reposition(2 * n_dofs, n_dofs, n_dofs, n_dofs);
    KIIR.reposition(2 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KII3R.reposition(2 * n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    KIuFR.reposition(2 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KIIFR.reposition(2 * n_dofs, 7 * n_dofs, n_dofs, n_dofs);
    //
    KI3uR.reposition(3 * n_dofs, 0, n_dofs, n_dofs);
    KI3nR.reposition(3 * n_dofs, n_dofs, n_dofs, n_dofs);
    KI3IR.reposition(3 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    KI3I3R.reposition(3 * n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    KI3uFR.reposition(3 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KI3I3FR.reposition(3 * n_dofs, 8 * n_dofs, n_dofs, n_dofs);
    //
    KCuR.reposition(4 * n_dofs, 0, n_dofs, n_dofs);
    KCCR.reposition(4 * n_dofs, 4 * n_dofs, n_dofs, n_dofs);
    KCuFR.reposition(4 * n_dofs, 5, n_dofs, n_dofs);
    KCCFR.reposition(4 * n_dofs, 9 * n_dofs, n_dofs, n_dofs);
    
    //
    KuuI.reposition(5 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KunI.reposition(5 * n_dofs, 6 * n_dofs, n_dofs, n_dofs);
    KuII.reposition(5 * n_dofs, 7 * n_dofs, n_dofs, n_dofs);
    KuI3I.reposition(5 * n_dofs, 8 * n_dofs, n_dofs, n_dofs);
    KuCI.reposition(5 * n_dofs, 9 * n_dofs, n_dofs, n_dofs);
    //
    KnuI.reposition(6 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KnnI.reposition(6 * n_dofs, 6 * n_dofs, n_dofs, n_dofs);
    KnII.reposition(6 * n_dofs, 7 * n_dofs, n_dofs, n_dofs);
    KnI3I.reposition(6 * n_dofs, 8 * n_dofs, n_dofs, n_dofs);
    KnuFI.reposition(6 * n_dofs, 0, n_dofs, n_dofs);
    KnnFI.reposition(6 * n_dofs, n_dofs, n_dofs, n_dofs);
    //
    KIuI.reposition(7 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KInI.reposition(7 * n_dofs, 6 * n_dofs, n_dofs, n_dofs);
    KIII.reposition(7 * n_dofs, 7 * n_dofs, n_dofs, n_dofs);
    KII3I.reposition(7 * n_dofs, 8 * n_dofs, n_dofs, n_dofs);
    KIuFI.reposition(7 * n_dofs, 0, n_dofs, n_dofs);
    KIIFI.reposition(7 * n_dofs, 2 * n_dofs, n_dofs, n_dofs);
    //
    KI3uI.reposition(8 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KI3nI.reposition(8 * n_dofs, 6 * n_dofs, n_dofs, n_dofs);
    KI3II.reposition(8 * n_dofs, 7 * n_dofs, n_dofs, n_dofs);
    KI3I3I.reposition(8 * n_dofs, 8 * n_dofs, n_dofs, n_dofs);
    KI3uFI.reposition(8 * n_dofs, 0, n_dofs, n_dofs);
    KI3I3FI.reposition(8 * n_dofs, 3 * n_dofs, n_dofs, n_dofs);
    //
    KCuI.reposition(9 * n_dofs, 5 * n_dofs, n_dofs, n_dofs);
    KCCI.reposition(9 * n_dofs, 9 * n_dofs, n_dofs, n_dofs);
    KCuFI.reposition(9 * n_dofs, 0, n_dofs, n_dofs);
    KCCFI.reposition(9 * n_dofs, 4 * n_dofs, n_dofs, n_dofs);
    //
    // Source Vector
    FuR.reposition(0, n_dofs);
    FnR.reposition(n_dofs, n_dofs);
    FIR.reposition(2 * n_dofs, n_dofs);
    FI3R.reposition(3 * n_dofs, n_dofs);
    FCR.reposition(4 * n_dofs, n_dofs);
    FuI.reposition(5 * n_dofs, n_dofs);
    FnI.reposition(6 * n_dofs, n_dofs);
    FII.reposition(7 * n_dofs, n_dofs);
    FI3I.reposition(8 * n_dofs, n_dofs);
    FCI.reposition(9 * n_dofs, n_dofs);
    //
    // Solution Vector
    XuR.reposition(0, n_dofs);
    XnR.reposition(n_dofs, n_dofs);
    XIR.reposition(2 * n_dofs, n_dofs);
    XI3R.reposition(3 * n_dofs, n_dofs);
    XCR.reposition(4 * n_dofs, n_dofs);
    XuI.reposition(5 * n_dofs, n_dofs);
    XnI.reposition(6 * n_dofs, n_dofs);
    XII.reposition(7 * n_dofs, n_dofs);
    XI3I.reposition(8 * n_dofs, n_dofs);
    XCI.reposition(9 * n_dofs, n_dofs);
    //
    // Solution Steady state Vector
    Su.reposition(0, n_dofs);
    Sn.reposition(n_dofs, n_dofs);
    SI.reposition(2 * n_dofs, n_dofs);
    SI3.reposition(3 * n_dofs, n_dofs);
    SC.reposition(4 * n_dofs, n_dofs);



    DSSCModel* sc =
      dynamic_cast<DSSCModel*>(
          device.get_material(subdomain)->get_model(get_id()));

    assert(sc != NULL);
    sc->reinit(elem);

    // Get the temperature given the element
    //vector<double> T_nodes = sc->get_temperature_at_nodes();

    vector<vector<double> > local_scaling(elem->n_nodes(), vector<double>(5, 1));
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
        u  += phi[i][qp] * Su(i);
        en += phi[i][qp] * Sn(i);
        eI += phi[i][qp] * SI(i);
        eI3 += phi[i][qp] * SI3(i);
        eC += phi[i][qp] * SC(i);
        e_field += dphi[i][qp] * Su(i);
        grad_en += dphi[i][qp] * Sn(i);
        grad_eI += dphi[i][qp] * SI(i);
        grad_eI3 += dphi[i][qp] * SI3(i);
        grad_eC += dphi[i][qp] * SC(i);
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
      sc->calculate_traps();
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
      //double sigma_n = (mu_n * n_e ) / C0_e;
      double sigma_n = (mu_n * n_e ) / C0_e;
      double sigma_I = (mu_I * n_I ) / C0_I;
      double sigma_I3 = (mu_I3 * n_I3 ) / C0_I3;
      double sigma_C = (mu_C * n_C ) / C0_C;

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
          double laplace = J * (dphi[i][qp] * dphi[j][qp]);

          KuuR(i,j) += l2_eps * laplace / local_scaling[i][4];
          KuuI(i,j) += l2_eps * laplace / local_scaling[i][4];


          if (!poisson_only())
          {
            if (is_TiO2)
            {
              KnnR(i,j) += sigma_n * laplace / local_scaling[i][0];
              KnnI(i,j) += sigma_n * laplace / local_scaling[i][0];
            }

            if (is_electrolyte)
            {
              KIIR(i,j) += sigma_I * laplace / local_scaling[i][1];
              KIII(i,j) += sigma_I * laplace / local_scaling[i][1];

              KI3I3R(i,j) += sigma_I3 * laplace / local_scaling[i][2];
              KI3I3I(i,j) += sigma_I3 * laplace / local_scaling[i][2];

              KCCR(i,j) += sigma_C * laplace / local_scaling[i][3];
              KCCI(i,j) += sigma_C * laplace / local_scaling[i][3];
            }
          }
        }

        if (!is_TiO2 || poisson_only())
        {
          KnnR(i,i) += 1;
          KnnI(i,i) += 1;
        }

        if (!is_electrolyte || poisson_only())
        {
          KIIR(i,i) += 1;
          KIII(i,i) += 1;
          KI3I3R(i,i) += 1;
          KI3I3I(i,i) += 1;
          KCCR(i,i) += 1;
          KCCI(i,i) += 1;
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
        double dn_dphi_exp_DOS = sc->get_density_derivative_n_exp_DOS();

        double drho = dC_dphi - (dn_dphi_exp_DOS + dI_dphi + dI3_dphi);
        drho *= phi0 / C0_tot;


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
        //double n_eq = sc->get_equilibrium_concentrations().n;
        //double dsigma_n = J * phi0 / C0_e * mu_n *  (1.4 * pow(n_e , 0.4)) * dn_dphi;
        double dsigma_n = J * phi0 / C0_e * mu_n * dn_dphi;
        double dsigma_I = J * phi0 / C0_I * mu_I * dI_dphi;
        double dsigma_I3 = J * phi0 / C0_I3 * mu_I3 * dI3_dphi;
        double dsigma_C = J * phi0 / C0_C * mu_C * dC_dphi;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          //double lap_e = (dphi[i][qp] * grad_en) / local_scaling[i][0];
          //double lap_I = (dphi[i][qp] * grad_eI) / local_scaling[i][1];
          //double lap_I3 = (dphi[i][qp] * grad_eI3) / local_scaling[i][2];
          //double lap_C = (dphi[i][qp] * grad_eC) / local_scaling[i][3];
          double lap_e = (dphi[i][qp] * grad_en);
          double lap_I = (dphi[i][qp] * grad_eI);
          double lap_I3 = (dphi[i][qp] * grad_eI3);
          double lap_C = (dphi[i][qp] * grad_eC);
          double dsigma_n_x_lap = dsigma_n * lap_e;
          double dsigma_I_x_lap = dsigma_I * lap_I;
          double dsigma_I3_x_lap = dsigma_I3 * lap_I3;
          double dsigma_C_x_lap = dsigma_C * lap_C;

          double volume = elem->volume();

          cons_cat->add(dof_indices_eC_R[i], -phi0 * J * dC_dphi * phi[i][qp] / C0_C / scaling_C );
          cons_cat->add(dof_indices_u_R[i], phi0 * J * dC_dphi * phi[i][qp]  / C0_C / scaling_C );

          cons_iodine->add(dof_indices_eI_R[i], -phi0 * J * dI_dphi * phi[i][qp] / (3.0*C0_tot*scaling_tot) );
          cons_iodine->add(dof_indices_eI3_R[i], -phi0 * J * dI3_dphi * phi[i][qp] / (C0_tot * scaling_tot) );
          cons_iodine->add(dof_indices_u_R[i], phi0 * J * phi[i][qp] * ( ( (1/3.0) * dI_dphi + dI3_dphi ) / (C0_tot * scaling_tot) ) );

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
              KnuR(i,j) += dsigma_n_x_phi / local_scaling[i][0];
              KnuI(i,j) += dsigma_n_x_phi / local_scaling[i][0];
              KnnR(i,j) -= dsigma_n_x_phi / local_scaling[i][0];
              KnnI(i,j) -= dsigma_n_x_phi / local_scaling[i][0];
            }

            if (is_electrolyte && !poisson_only())
            {
              KIuR(i,j) += dsigma_I_x_phi / local_scaling[i][1];
              KIuI(i,j) += dsigma_I_x_phi / local_scaling[i][1];
              KIIR(i,j) -= dsigma_I_x_phi / local_scaling[i][1];
              KIII(i,j) -= dsigma_I_x_phi / local_scaling[i][1];

              KI3uR(i,j) += dsigma_I3_x_phi / local_scaling[i][2];
              KI3uI(i,j) += dsigma_I3_x_phi / local_scaling[i][2];
              KI3I3R(i,j) -= dsigma_I3_x_phi / local_scaling[i][2];
              KI3I3I(i,j) -= dsigma_I3_x_phi / local_scaling[i][2];

              KCuR(i,j) += dsigma_C_x_phi / local_scaling[i][3];
              KCuI(i,j) += dsigma_C_x_phi / local_scaling[i][3];
              KCCR(i,j) -= dsigma_C_x_phi / local_scaling[i][3];
              KCCI(i,j) -= dsigma_C_x_phi / local_scaling[i][3];
            }




            // The dFe_i/dX_j part
            double phi_i_x_phi_j = J * phi[i][qp] * phi[j][qp];

            KuuR(i,j) -= drho * phi_i_x_phi_j / local_scaling[i][4];
            KuuI(i,j) -= drho * phi_i_x_phi_j / local_scaling[i][4];

            if (is_TiO2)
            {
              KunR(i,j) -= phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KunI(i,j) -= phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;

              if (!poisson_only())
              {
                KnuR(i,j) -= dR[0] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                KnuI(i,j) -= dR[0] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                KnnR(i,j) -= dR[1] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                KnnI(i,j) -= dR[1] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
              
                KnuFR(i,j) += -freq * phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][0] / C0_e;
                KnnFR(i,j) +=  freq * phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][0] / C0_e;
                
                KnuFI(i,j) +=  freq * phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][0] / C0_e;
                KnnFI(i,j) += -freq * phi0 * dn_dphi_exp_DOS * phi_i_x_phi_j / local_scaling[i][0] / C0_e;
              }
            }

            if (is_electrolyte)
            {
              KuIR(i,j) -= phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuII(i,j) -= phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuI3R(i,j) -= phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuI3I(i,j) -= phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuCR(i,j) += phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;
              KuCI(i,j) += phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][4] / C0_tot;

            
              //Kau(i,j) -= phi0 * dC_dphi * phi_i_x_phi_j / C0_C;
              //KaC(i,j) += phi0 * dC_dphi * phi_i_x_phi_j / C0_C;

              //Kbu(i,j) -= phi0 * (dI3_dphi + dI_dphi / 3.0) * phi_i_x_phi_j / C0_tot;
              //KbI(i,j) += phi0 * dI_dphi / 3.0 * phi_i_x_phi_j / C0_tot;
              //KbI3(i,j) += phi0 * dI3_dphi * phi_i_x_phi_j / C0_tot;

              if (!poisson_only())
              {
                KIuFR(i,j) += -freq * phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][1] / C0_I;
                KIIFR(i,j) +=  freq * phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][1] / C0_I;
                
                KIuFI(i,j) +=  freq * phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][1] / C0_I;
                KIIFI(i,j) += -freq * phi0 * dI_dphi * phi_i_x_phi_j / local_scaling[i][1] / C0_I;
            
                KI3uFR(i,j) += -freq * phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I3;
                KI3I3FR(i,j) +=  freq * phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I3;
                
                KI3uFI(i,j) +=  freq * phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I3;
                KI3I3FI(i,j) += -freq * phi0 * dI3_dphi * phi_i_x_phi_j / local_scaling[i][2] / C0_I3;
            
                KCuFR(i,j) += -freq * phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][3] / C0_C;
                KCCFR(i,j) +=  freq * phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][3] / C0_C;
                
                KCuFI(i,j) +=  freq * phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][3] / C0_C;
                KCCFI(i,j) += -freq * phi0 * dC_dphi * phi_i_x_phi_j / local_scaling[i][3] / C0_C;

                if (is_TiO2)
                {
                  KnIR(i,j) -= dR[2] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                  KnII(i,j) -= dR[2] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                  KnI3R(i,j) -= dR[3] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;
                  KnI3I(i,j) -= dR[3] * phi_i_x_phi_j / local_scaling[i][0] / R0_e;

                  KInR(i,j) += 1.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                  KInI(i,j) += 1.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                  KI3nR(i,j) -= 0.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                  KI3nI(i,j) -= 0.5 * dR[1] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                }

                KIuR(i,j) += 1.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KIuI(i,j) += 1.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KIIR(i,j) += 1.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KIII(i,j) += 1.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII3R(i,j) += 1.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;
                KII3I(i,j) += 1.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][1] / R0_I;

                KI3uR(i,j) -= 0.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3uI(i,j) -= 0.5 * dR[0] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3IR(i,j) -= 0.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3II(i,j) -= 0.5 * dR[2] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3I3R(i,j) -= 0.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
                KI3I3I(i,j) -= 0.5 * dR[3] * phi_i_x_phi_j / local_scaling[i][2] / R0_I3;
              }
            }
          }

        }

      }


      // if we are doing residual, calculate rhs contribution (i.e. Fe)
      if (residual != NULL)
      {
        // charge density
        double J_x_rho = J * sc->get_charge_density() / C0_tot;
        //double J_x_rho = 0.0;

        // net recombination rate
        double J_x_R = 0.0;
        if ((is_TiO2 && is_electrolyte) && !poisson_only())
          J_x_R = J * R;


        for (unsigned int i = 0; i < n_dofs; i++)
        {
          double net_recomb = J_x_R * phi[i][qp];

//          Fu(i) -= J_x_rho * phi[i][qp] / local_scaling[i][4];

//          Fn(i) -= net_recomb / R0_e / local_scaling[i][0];

          // TODO factors
//          FI(i) += 1.5 * net_recomb / R0_I / local_scaling[i][1];
//          FI3(i) -= 0.5 * net_recomb / R0_I3 / local_scaling[i][2];

          double volume = elem->volume();

          tot_cat += J * phi[i][qp] * n_C / C0_C;
          tot_iodine += J * phi[i][qp] * (n_I3 + n_I / 3.0) / C0_tot;

        }
      }

    } // end loop over quadrature points


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
        // we need to know if it is an outer or inner boundary
        bool true_boundary = environment.is_outer_boundary(side);
        bool inner_boundary = environment.is_inner_boundary(side);

        Boundary* boundary = environment.get_boundary(side);

        DSSCContact* contact = NULL;
        if (boundary != NULL)
          contact = dynamic_cast<DSSCContact*>(
              boundary->get_boundary_properties(get_id()));

        AutoPtr<Elem> side(elem->build_side(s));

        fe_face->reinit(elem, s);



        // calculate the fluxes on the nodes
        if ((contact != NULL) && (dim > 1))
        {
          vector<Point> p(side->n_nodes());
          for (unsigned int i = 0; i < side->n_nodes(); i++)
          {
            nodes_on_boundary_sides.insert(side->node(i));
            p[i] = side->point(i);
          }
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

         if (true_boundary)
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
                u  += phi_face[i][qp] * Su(i);
                en += phi_face[i][qp] * Sn(i);
                eI += phi_face[i][qp] * SI(i);
                eI3 += phi_face[i][qp] * SI3(i);
                eC += phi_face[i][qp] * SC(i);
                e_field += dphi_face[i][qp] * Su(i);
                grad_en += dphi_face[i][qp] * Sn(i);
                grad_eI += dphi_face[i][qp] * SI(i);
                grad_eI3 += dphi_face[i][qp] * SI3(i);
                grad_eC += dphi_face[i][qp] * SC(i);
              }

              sc->set_potentials(phi0 * u, phi0 * en, phi0 * eI, phi0 * eI3, phi0 * eC);
              sc->set_coordinates(q_point_face[qp]);

              sc->calculate_densities();
              sc->calculate_traps();
              sc->calculate_equilibrium_traps();

              // we put the phi0 here for convenience
              double sigma_I = sc->get_mobility_I() * sc->get_density_I() / C0_I;
              double dsigma_I = phi0 / C0_I * sc->get_mobility_I() *
                sc->get_density_derivative_I();
              double sigma_I3 = sc->get_mobility_I3() * sc->get_density_I3() / C0_I3;
              double dsigma_I3 = phi0 / C0_I3 * sc->get_mobility_I3() *
                sc->get_density_derivative_I3();
              double sigma_n = sc->get_mobility_n() * sc->get_density_n() / C0_e;
              double dsigma_n = phi0 / C0_e * sc->get_mobility_n() *
                sc->get_density_derivative_n();
              double sigma_C = sc->get_mobility_C() * sc->get_density_C() / C0_C;
              double dsigma_C = phi0 / C0_C * sc->get_mobility_C() *
                sc->get_density_derivative_C();


            //define flux direction
            

            // the jacobian x weight x scaling
            double J = JxW_face[qp];

            if (contact->is_cathode())
            {
              for (unsigned int i = 0; i < n_dofs; i++)
              {

                if (jacobian != NULL)
                {
                  for (unsigned int j = 0; j < n_dofs; j++)
                  {

                    //double res = contact->get_load() * x0;
                    double j0 = contact->get_ex_curr();
                
                    double kT = sc->get_lattice_temperature();
                    double kT2 = 2*kT;
	            //double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	            //double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
                    double Normal_I = x0 / (phi0 * C0_I * Constants::e * scaling_I);
                    double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * scaling_I3);
                
                    //double Eredox = -2*u + 1.5*eI - 0.5*eI3;
                    double Eredox = -u;

                    double J_phi_i_phi_j = J * phi_face[i][qp] * phi_face[j][qp];
	        
                    //KIu(i,j) += 1.5 * 2 * J_phi_i_phi_j * Normal_I / res;
                    //KII(i,j) += 1.5 * -1.5 * J_phi_i_phi_j * Normal_I / res;
                    //KII3(i,j) += 1.5 * 0.5 * J_phi_i_phi_j * Normal_I / res;
//                    KIu(i,j) += 1.5 * J_phi_i_phi_j * Normal_I / res;
                    KIuR(i,j) += 1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I / kT2;
                    KIuI(i,j) += 1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I / kT2;

                    //KI3u(i,j) += -0.5 * 2 * J_phi_i_phi_j * Normal_I3 / res;
                    //KI3I(i,j) += -0.5 * -1.5 * J_phi_i_phi_j * Normal_I3 / res;
                    //KI3I3(i,j) += -0.5 * 0.5 * J_phi_i_phi_j * Normal_I3 / res;
//                    KI3u(i,j) += -0.5 * J_phi_i_phi_j * Normal_I3 / res;
                    KI3uR(i,j) += -0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I3 / kT2;
                    KI3uI(i,j) += -0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * J_phi_i_phi_j * Normal_I3 / kT2;

                  }
                }
                if (residual != NULL)
                {
                    
                }
              }

            }
            else
            {

             if (n_cat == NULL)
             {
               n_cat = side->get_node(0);
             }

            for (unsigned int i = 0; i < n_dofs; i++)
            {
               if (jacobian != NULL)
               {
               for (unsigned int j = 0; j < n_dofs; j++)
                {
	           //double Normal_n = x0 / (phi0 * C0_e * local_scaling[s][0] );
	           double Normal_n = x0 / (phi0 * C0_e );
                   double kinetic_anode = contact->get_kinetic();
                    
                   double J_phi_i_phi_j = J * phi_face[i][qp] * phi_face[j][qp];

                   double dn_dphi = sc->get_density_derivative_n();

                   KnuR(i,j) += -kinetic_anode * dn_dphi * Normal_n * J_phi_i_phi_j;
                   KnuI(i,j) += -kinetic_anode * dn_dphi * Normal_n * J_phi_i_phi_j;
                   KnnR(i,j) += -kinetic_anode * -dn_dphi * Normal_n * J_phi_i_phi_j;
                   KnnI(i,j) += -kinetic_anode * -dn_dphi * Normal_n * J_phi_i_phi_j;
                }
               }
               if (residual != NULL)
               { 
                 
                   double density_n = sc->get_density_n();
                   double n_dark = sc->get_equilibrium_concentrations().n;
                   double bias = contact->get_potential();
                   double kT = sc->get_lattice_temperature();
	           //double Normal_n = x0 / ( phi0 * C0_e * Constants::e * local_scaling[s][0] );
	           double Normal_n = x0 / ( phi0 * C0_e * local_scaling[s][0] );
                   double kinetic_anode = contact->get_kinetic();

                   FnR(s) += -kinetic_anode * n_dark * Normal_n / kT; 
               }
             }

            }

         } //contact

        } // qp

        } // true_boundary
        } // dimension
        else // i.e. dim == 1
        {


        if (true_boundary)
        {
          nodes_on_boundary_sides.insert(elem->node(s));

          // s is the node of the element lying on the boundary
          Real u  = Su(s);
          Real en = Sn(s);
          Real eI = SI(s);
          Real eI3 = SI3(s);
          Real eC = SC(s);

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
            grad_en += dphi_face[n][0](0) * Sn(n);
            grad_eI += dphi_face[n][0](0) * SI(n);
            grad_eI3 += dphi_face[n][0](0) * SI3(n);
            grad_eC += dphi_face[n][0](0) * SC(n);
            e_field(0) += dphi_face[n][0](0) * Su(n);
          }

          sc->set_electric_field(phi0 / x0 * e_field);
          sc->set_grad_fermi_n(phi0 / x0 * RealGradient(grad_en, 0.0, 0.0));
          sc->set_grad_fermi_I(phi0 / x0 * RealGradient(grad_eI, 0.0, 0.0));
          sc->set_grad_fermi_I3(phi0 / x0 * RealGradient(grad_eI3, 0.0, 0.0));
          sc->set_grad_fermi_C(phi0 / x0 * RealGradient(grad_eC, 0.0, 0.0));
          sc->calculate_densities();
          sc->calculate_traps();
          sc->calculate_equilibrium_traps();

          // we put the phi0 here for convenience
          double sigma_I = sc->get_mobility_I() * sc->get_density_I() / C0_I;
          double dsigma_I = phi0 / C0_I * sc->get_mobility_I() *
            sc->get_density_derivative_I();
          double sigma_I3 = sc->get_mobility_I3() * sc->get_density_I3() / C0_I3;
          double dsigma_I3 = phi0 / C0_I3 * sc->get_mobility_I3() *
            sc->get_density_derivative_I3();
          double sigma_n = sc->get_mobility_n() * sc->get_density_n() / C0_e;
          double dsigma_n = phi0 / C0_e * sc->get_mobility_n() *
            sc->get_density_derivative_n();
          double sigma_C = sc->get_mobility_C() * sc->get_density_C() / C0_C;
          double dsigma_C = phi0 / C0_C * sc->get_mobility_C() *
            sc->get_density_derivative_C();

          double x_c = elem->centroid()(0);
          double x_s = elem->point(s)(0);
          double sign = (x_s > x_c) ? 1 : -1;


          if (contact->is_cathode())
          {
//             if (n_cat == NULL)
//             {
//               n_cat = side->get_node(0);            
//             }

            if (jacobian != NULL)
             {

                //double res = contact->get_load() * x0;
                
                double j0 = contact->get_ex_curr();
                double kT = sc->get_lattice_temperature();
                double kT2 = 2*kT;
	        double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	        double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
                
                double Eredox = -Su(s);
               
                double bias = contact->get_potential();
                
                double Eredox1 = bias + 2.5*SI(s) - 1.5*SI3(s) - 2*Su(s);
                double Eredox2 = -bias +4*Su(s) - 3.5*SI(s) + 0.5*SI3(s);

//                KI3u(s,s) += 0.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//                KI3I(s,s) += 0.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//                KI3I3(s,s) += 0.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
                
//                KIu(s,s) += -1.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//                KII(s,s) += -1.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//                KII3(s,s) += -1.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
	        

//               KIu(s,s) += 1.5 * Normal_I / res;
//               KI3u(s,s) += -0.5 * Normal_I3 / res;
                
                KIuR(s,s) += -1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I / kT2;
                KIuI(s,s) += -1.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I / kT2;

                KI3uR(s,s) += 0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I3 / kT2;
                KI3uI(s,s) += 0.5 * j0 * ( exp(Eredox/kT2) + exp(-Eredox/kT2) ) * (-1) * Normal_I3 / kT2;
                

            }
            if (residual != NULL)
            {


             }
           }
           else
	   {

             if (n_cat == NULL)
             {
               n_cat = side->get_node(0);            
             }

             if (jacobian != NULL)
             {
               
	       //double Normal_n = x0 / (phi0 * C0_e * Constants::e * local_scaling[s][0] );
	       double Normal_n = x0 / (phi0 * C0_e * local_scaling[s][0] );
               double kinetic_anode = contact->get_kinetic();
               
               double bias = contact->get_potential();
               double kT = sc->get_lattice_temperature();
               double j0 = contact->get_ex_curr();
               double kT2 = 2*kT;
	       double Normal_I = x0 / (phi0 * C0_I * Constants::e * local_scaling[s][1] );
	       double Normal_I3 = x0 / (phi0 * C0_I3 * Constants::e * local_scaling[s][2] );
               
               double Eredox1 = bias + 2.5*SI(s) - 1.5*SI3(s) - 2*Su(s);
               double Eredox2 = -bias +4*Su(s) - 3.5*SI(s) + 0.5*SI3(s);
                 
               double dn_dphi = sc->get_density_derivative_n();
               
//               KI3u(s,s) += 0.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//               KI3I(s,s) += 0.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//               KI3I3(s,s) += 0.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;
               
//               KIu(s,s) += -1.5 * j0 * ( -2*exp( (Eredox1)/kT2 ) -4*exp( (Eredox2)/kT2 ) ) * Normal_I3/kT2;
//               KII(s,s) += -1.5 * j0 * ( 2.5*exp( (Eredox1)/kT2 ) + 3.5*exp( (Eredox2)/kT2 ) )  * Normal_I3 / kT2;
//               KII3(s,s) += -1.5 * j0 * ( -1.5*exp( (Eredox1)/kT2 ) -0.5*exp( -(Eredox2)/kT2 ) ) * Normal_I3 /kT2;

               KnuR(s,s) += -kinetic_anode * dn_dphi * Normal_n;
               KnuI(s,s) += -kinetic_anode * dn_dphi * Normal_n;
               KnnR(s,s) += -kinetic_anode * -dn_dphi * Normal_n;
               KnnI(s,s) += -kinetic_anode * -dn_dphi * Normal_n;

             }
             if (residual != NULL)
             { 
                 
               double density_n = sc->get_density_n();
               double n_dark = sc->get_equilibrium_concentrations().n;
               double bias = contact->get_potential();
               double kT = sc->get_lattice_temperature();
               double n_0 = n_dark * exp(bias/kT); 
	       //double Normal_n = x0 / ( phi0 * C0_e * Constants::e * local_scaling[s][0] );
	       double Normal_n = x0 / ( phi0 * C0_e * local_scaling[s][0] );
               double kinetic_anode = contact->get_kinetic();

               FnR(s) += -kinetic_anode * n_dark * Normal_n / kT; 

             }

           }
	}
       } // if (dim ...)
      }
    } // end loop over element sides


    // check if it is a dielectric
    if (!sc->is_TiO2())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          KunR(i, j) = KunI(i, j) = 0.0;
          KnuR(i, j) = KnuI(i, j) = KnnR(i, j) = KnnI(i, j) = KnIR(i, j) = KnII(i, j) = KnI3R(i, j) = KnI3I(i, j) = 0.0;
          KInR(i, j) = KInI(i, j) = KI3nR(i, j) = KI3nI(i, j) = 0.0;
        }

        if (!is_internal_boundary_node(elem->get_node(i)))
          KnnR(i, i) = KnnI(i, i) = 1.0;

        // we simply set the electrochemical potential to zero
        FnR(i) = FnI(i) = 0.0;
      }
    }

    // check if it is a dielectric
    if (!sc->is_electrolyte())
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        for (unsigned int j = 0; j < n_dofs; j++)
        {
          KuII(i, j) = KuIR(i, j) = KuI3R(i, j) = 0.0; 
          KuI3I(i, j) = KuCR(i, j) = KuCI(i, j) = 0.0;

          KIuR(i, j) = KIuI(i, j) = KInR(i, j) = KInI(i, j) = 0.0;
          KIIR(i, j) = KIII(i, j) = KII3R(i, j) = KII3I(i, j) = 0.0;
          KI3uR(i, j) = KI3uI(i, j) = KI3nR(i, j) = KI3nI(i, j) = 0.0;
          KI3IR(i, j) = KI3II(i, j) = KI3I3R(i, j) = KI3I3I(i, j) = 0.0;

          KCuR(i, j) = KCuI(i, j) = KCCR(i, j) = KCCI(i, j) = 0.0;
          
          KnIR(i, j) = KnII(i, j) = KnI3R(i, j) = KnI3I(i, j) = 0.0;
        }

        if (!is_internal_boundary_node(elem->get_node(i)))
        {
          KIIR(i, i) = KIII(i, i) = 1.0;
          KI3I3R(i, i) = KI3I3I(i, i) = 1.0;
          KCCR(i, i) = KCCI(i, i) = 1.0;
        }

        // we simply set the electrochemical potential to zero
        FIR(i) = FII(i) = 0.0;
        FI3R(i) = FI3I(i) = 0.0;
        FCR(i) = FCI(i) = 0.0;
      }
    } 

    // constrain the jacobian and the rhs to account for constrained
    // DOFs
    // NOTE: this changes dof_indices that's why th application of
    //       Dirichlet type BCs needs special care
    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

///*
    //
    // now as last thing we apply Dirichlet type Bcs
    //
    BoundaryNodeList::const_iterator node_it;
    const BoundaryNodeList::const_iterator end =
      dirichlet_nodes.end();
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

          // we only impose Dirichlet type BCs if the node has an associated
          // boundary side
          if (nodes_on_boundary_sides.find(elem->node(i)) !=
              nodes_on_boundary_sides.end())
          {
            if ( (!contact->is_cathode()) && (!contact->is_gate()) )
            {

              double bias = contact->get_potential();

//              Ke.condense(i + n_dofs, i + n_dofs, bias, Fe);

            }
            else if (contact->is_gate())
            {
               double bias = contact->get_potential();

              Ke.condense(i, i, bias, Fe);
            }
          }
        }
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

  //COMMENT: dovrei mettere la conservazione due volte

//  unsigned int dof_cat = n_cat->dof_number(system.number(), eC_var, 0);
//  unsigned int dof_iodine = n_cat->dof_number(system.number(), eI_var, 0);
  unsigned int dof_cat_R = n_cat->dof_number(system.number(), eC_var_R, 0);
  unsigned int dof_iodine_R = n_cat->dof_number(system.number(), eI_var_R, 0);
  unsigned int dof_cat_I = n_cat->dof_number(system.number(), eC_var_I, 0);
  unsigned int dof_iodine_I = n_cat->dof_number(system.number(), eI_var_I, 0);

  if (jacobian != NULL)
  {

    jacobian->close();

    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_cat->size(); j++)
      jacobian->set(dof_cat_R, j, (*cons_cat)(j));

    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_iodine->size(); j++)
      jacobian->set(dof_iodine_R, j, (*cons_iodine)(j));

    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_cat->size(); j++)
      jacobian->set(dof_cat_I, j, (*cons_cat)(j));
    
    assert(n_cat != NULL);
    for (unsigned int j = 0; j < cons_cat->size(); j++)
      jacobian->set(dof_cat_I, j, (*cons_cat)(j));
    
    jacobian->close();
    //jacobian->print_matlab("J.m");
  }
  else
  {
    residual->close();
    residual->set(dof_cat_R, (tot_cat / scaling_C - _cation_amount / C0_C / scaling_C));
    residual->set(dof_iodine_R, (tot_iodine / scaling_tot - _iodine_amount / C0_tot / scaling_tot));
    residual->close();
    //residual->print_matlab("F.m");
  }


  perf_log.stop_event("assembly frequency");

}





