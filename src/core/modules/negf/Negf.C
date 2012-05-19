// $Id: Negf.h 2964 2011-10-10 20:34:57Z fpalomba $

#include "Negf.h"
#include "NegfModel.h"
#include "Device.h"
#include "QuantumContact.h"
#include "SimulationEnvironment.h"
#include "SimulationOptions.h"
#include "DataOutput.h"
#include "Messages.h"
#include "InitFailedException.h"
#include "PotentialInterface.h"
#include "KspaceIntegration.h"
#include "KspaceIntegrationTemplate.h"

// Basic include files needed for the mesh functionality.
#include "fe.h"
#include "fe_interface.h"
// Define generic quadrature rules.
#include "quadrature.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "mesh.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"
#include "tensor_value.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_macro.h"
#include "petscversion.h"


// C++ includes
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cassert>
#include <math.h>

using namespace Constants;

TIBER_MODULE(Negf, MODULE_NAME)

Negf* Negf::static_this;

Negf::Negf(const ModelOptions& options) :
                              SimulationInterface(options)
{
  _device_n_dofs = 0;
  _k_int_density = NULL;
  _k_int_current = NULL;
  _libnegf = NegfWrapper::create();

  this->has_solution_vector(false);
}

Negf::~Negf(void)
{
  delete _libnegf;
}

Negf*
Negf::create(const ModelOptions& options)
{
  static_this = new Negf(options);
  std::cout<<"this & "<<static_this<<std::endl;
  return static_this;
}

PhysicalModel*
Negf::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  return  NegfModel::create(mat, options);

}


//! The initialization
void
Negf::do_init(void)
{
  static_this = this;

  parse_options();

  _env = &get_environment();

  _device = &(_env->get_device());

  // Prepare QuantumContact Map
  // Note: Boundary _contact_names are the same as QC names !!
  {
    _quantum_contacts.clear();
    SimulationEnvironment::BoundaryIterator it = _env->boundaries_begin();
    const SimulationEnvironment::BoundaryIterator end = _env->boundaries_end();
    for (; it != end; ++it)
    {
      QuantumContact* qc = _device->get_quantum_contact((*it)->get_name());

      if (qc != NULL)
      {
        _quantum_contacts[qc->get_id()] = qc;
        _qc_boundaries[*it] = qc;
        _bd_map[qc] = *it;

        std::cerr<<"_quantum_contact "<<qc->get_id()<<" "<<(*it)->get_name()<<std::endl;
        // Quantum Contacts are activated here: so dof_map come out correctly
        qc->activate_elements();
        qc->set_neighbor_map();

      }
      else
      {
        if ((*it)->get_options().get_option("dirichlet",false))
        {
          _dirichlet_boundaries.insert(*it);
          std::cout<<"dirichlet bc on "<<(*it)->get_name()<<std::endl;
        }
      }

    }
  }


  std::cout<<"this & "<<static_this<<std::endl;

  // get the number of subbands.
   const MeshBase& mesh = get_mesh();
   MeshBase::const_element_iterator el = mesh.active_elements_begin();
   const Elem* elem = *el;

  NegfModel* negfmod = get_bulk_model<NegfModel>(elem);
  // Setup a simple effective mass Hamiltonian
  std::cout<<"(negf) create eq sys for H"<<std::endl;
  ID id = create_equation_system("linear","");
  _sys_H = &get_equation_system<TiberLinearSystem>(id);

  // attach a variable for each subband
  for (unsigned int band=0; band < negfmod->get_n_bands(); band++)
  {
    std::stringstream out;
    out<<"phi"<<band;
    _sys_H->add_variable(out.str(), FIRST, LAGRANGE);
  }

  _sys_H->add_matrix("Hi");

  // Setup Overlap matrix
  std::cout<<"(negf) create eq sys for S"<<std::endl;
  id = create_equation_system("linear","");
  _sys_S = &get_equation_system<TiberLinearSystem>(id);

  // attach a variable for each subband
  for (unsigned int band=0; band < negfmod->get_n_bands(); band++)
  {
    std::stringstream out;
    out<<"phi"<<band;
    _sys_S->add_variable(out.str(), FIRST, LAGRANGE);
  }

  _sys_S->add_matrix("Si");


  // attach system for reorder dofs
  std::cout<<"(negf) get eq sys for reorder"<<std::endl;
  id = create_equation_system("linear","");
  _sys = &get_equation_system<TiberLinearSystem>(id);

  // We add a second system just to contain the density
  if (plot_solution(elDensity))
  {
    std::cout<<"(negf) create eq sys for elDensity"<<std::endl;
    id = create_equation_system("linear","");
    _qdens_sys = &get_equation_system<TiberLinearSystem>(id);
    _qdens_sys->add_variable("edens", libMeshEnums::FIRST);
    //_qdens_sys->add_variable("hdens", libMeshEnums::FIRST);

  }
  std::cout<<"(negf) init H and S"<<std::endl;
  _sys_H->attach_assemble_function(ham_assemble);

  _sys_H->init();
  _sys_S->init();

  std::cout<<"(negf) init k-integration"<<std::endl;
  init_k_space_integration();

  std::cout<<"(negf) activate quantum contacts"<<std::endl;
  activate_quantum_contacts();

  std::cout<<"(negf) reorder dofs"<<std::endl;
  reorder(); // dof indices reorder

  std::cout<<"(negf) init done"<<std::endl;
}

void
Negf::do_reinit(void)
{
  static_this = this;
  //std::cout<<"(negf) clean up libnegf"<<std::endl;
   _libnegf->clean_libnegf();

   //std::cout<<"(negf) clear systems"<<std::endl;

}

void
Negf::init_k_space_integration(void)
{
  if (get_options().has_submodel("k_integration_density"))
  {
    //-----------------DENSITY ---------------------------------------------
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_density"));
    ModelOptions& kopts = it->second;

    kopts.set_option("mesh_units",get_mesh_units());
    kopts.set_option("k_space_dimension",3 - get_mesh().mesh_dimension());
    kopts.set_option("verbose", SimulationOptions::verbose() );


    _k_int_density = KspaceIntegrationTemplate<Negf>::create(this,kopts);

    if (_k_int_density == NULL)
      throw InitFailedException("Could not create k-integration");

    _k_int_density->init();
  }
  //------------------CURRENT --------------------------------------------------
  if (get_options().has_submodel("k_integration_current"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_current"));
    ModelOptions& kopts = it->second;

    kopts.set_option("mesh_units",get_mesh_units());
    kopts.set_option("k_space_dimension",3 - get_mesh().mesh_dimension());
    kopts.set_option("verbose", SimulationOptions::verbose() );

    _k_int_current = KspaceIntegrationTemplate<Negf>::create(this,kopts);

    if (_k_int_current == NULL)
      throw InitFailedException("Could not create k-integration");

    _k_int_current->init();

  }
}

void
Negf::setup_effectivemass_hamil()
{
  do_reinit();

  _sys_H->assemble();

  print_ham("matlab");

  print_Lib();

  _libnegf->init();

  if (opt.verbosity > 60) _libnegf->partition_info();

  _libnegf->set_verbose(opt.verbosity);

  _libnegf->set_write_ldos(opt.writeLDOS);

  _libnegf->set_write_tunn(true);

  _libnegf->set_iteration(1);

  //Messages::info("setting scratch path to "+SimulationOptions::scratch_path);
  _libnegf->set_scratch_path(SimulationOptions::scratch_path);

  _libnegf->set_output_path(get_output_directory());

  // set not to compute Device-Contact blocks
  _libnegf->device_contact_dm(0);

  // set reference contact at maximum electrochem pot.
  _libnegf->set_reference(1);

}

void
Negf::compute_current(void)
{
  _libnegf->set_verbose(opt.verbosity);

  //std::cout << _libnegf->current() << std::endl;

  current.clear();
  current.resize(2,0.0);

  //TODO: elCurrent. hlCurrent

  current[0] = _libnegf->current();
  current[1] = current[0];
}

void
Negf::do_solve(void)
{
  static_this = this;

  activate_quantum_contacts();

  //reorder(); // dof indices reorder

  if ( plot_solution("elDensity") )
  {
    Messages::info("Computing Density");

    _qdens_sys->init();


    if (get_options().has_submodel("k_integration_density"))
    {
      _which_integration = 1;

      _k_int_density->solve();

      density = _k_int_density->get_solution();

    }
    else
    {
      _k_vec.zero();

      setup_effectivemass_hamil();

      calculate_density("el");
    }

    Messages::info("Density done");

  }


  if ( plot_solution("Current") )
  {
    Messages::info("Computing Current");

    if (get_options().has_submodel("k_integration_current"))
    {
      _which_integration = 3;

      _k_int_current->solve();

      current = _k_int_current->get_solution();

    }
    else
    {
       _k_vec.zero();

       setup_effectivemass_hamil();

       compute_current();
    }

    double u = get_mesh_units();
    unsigned int dim = get_mesh().mesh_dimension();
    double area_factor;

    switch (dim)
    {
      case 1:
        area_factor = 1e14;  // k^2 is in 1/nm^2 => 1e14 1/cm^2
        break;
      case 2:
        area_factor = 1e7;
        break;
      case 3:
        area_factor = 1.0;
    }

    ID id = 0;
    double sign;
    std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
    const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
    for (; it != end; ++it)
    {
      if(_contact_potential[it->second] == mumax)
        sign = 1.0;
      else
        sign =-1.0;

      _contact_current[it->second] = sign * current[id] * area_factor;
      id++;
    }

    Messages::info("Current done");
    plot_globaldata();

  }

  deactivate_quantum_contacts();
}

void
Negf::calculate_for_k_point(const Point& k_point,
                                   DofField& field,
                                   double& error)
{
   for(short i=0;i<3;i++) _k_vec(i) = k_point(i);

   setup_effectivemass_hamil();

   if (_which_integration == 1)
   {
       calculate_density("el");
       NumericVector<Number>& qdens = *_qdens_sys->solution;

       field.clear();
       field.reserve(_device_n_dofs);

       for (unsigned int el=0; el < _device_n_dofs; el++)
         field.push_back(qdens(el));

       error = 0.0;

       return;
   }

   if (_which_integration == 2)
   {
       return;
   }


   if (_which_integration == 3)
   {
       field.clear();
       compute_current();
       //double curr = _libnegf->current();
       //std::cout<<"current: "<< curr<<std::endl;
       //_contact_potential[]
       //field.push_back(curr);
       //field.push_back(-curr);
       field = current;

       error = 0.0;
       return;
   }

}


void
Negf::calculate_density(const std::string& particle)
{
  double u = get_mesh_units();
  unsigned int dim = get_mesh().mesh_dimension();

  deactivate_quantum_contacts();

  // compute total number of n_vars n_dofs
  unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_dofs = _qc_n_dofs[_quantum_contacts.size()-1];

  int particle_id;

  if (particle == "el") particle_id = _qdens_sys->variable_number("edens");
  if (particle == "hl") particle_id = _qdens_sys->variable_number("hdens");

  std::vector<double> density(_device_n_dofs * n_vars, 0.0);

  //-------------------------------------------------------
  // HERE WE COMPUTE DENSITY!
  //-------------------------------------------------------
  _libnegf->density(density);

  // ---- debug ---------------------------------------------
  std::string out_file = "density.dat";
  std::fstream ff(out_file.c_str(),std::fstream::out);
  for (unsigned int i = 0; i < density.size(); i++)
        ff<< density[i] <<std::endl;
  ff.close();
  // ---------------------------------------------------------


  double equ;
  switch (dim)
  {
    case 1:
      equ = 1.0e14/(u*1.0e2); //  k^2 is in 1/nm^2 => 1e14 1/cm^2  / mesh_units -> cm
      break;
    case 2:
      equ = 1.0e7/(u*u*1.0e4); //
      break;
    case 3:
      equ = 1.0/(u*u*u*1.0e6);
  }

  DofMap& dof_map_qdens = _qdens_sys->get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  // setup output vector qdens
  NumericVector<Number>& qdens = *_qdens_sys->solution;
  qdens.zero();

  // compute connectivity of each node
  // we need the connectivity of the nodes to not double count
  //std::vector<int> connectivity(_device_n_dofs * n_vars, 0);
  //std::cout<<"qdens.size: "<<qdens.size()<<std::endl;

  std::vector<int> connectivity(qdens.size(), 0);
  {
    MeshBase::const_element_iterator el = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, particle_id);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        connectivity[dof_indices_qdens[n]]++;
    }
  }


  // The _sys_H system contains the nodal dofmap (all variables)
  DofMap& dof_map = _sys_H->get_dof_map();
  std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator el = get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    NegfModel* negfmod = get_bulk_model<NegfModel>(elem);

    dof_map_qdens.dof_indices(elem, dof_indices_qdens, particle_id);

    for (int band = 0; band < n_vars; band++)
    {
      dof_map.dof_indices(elem, dof_indices, band);

      double deg = negfmod->get_degeneracy(band);


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        double dens = density[_inv_perm[dof_indices[n]]];

        //std::cout<<band*n_dofs+n<<" dens: "<<dens<<std::endl;

        double val = equ * deg * dens / connectivity[dof_indices_qdens[n]];

        qdens.add( dof_indices_qdens[n], val);
      }

      //for (unsigned int n = 0; n < elem->n_nodes(); n++)
      //std::cout<<connectivity[ dof_indices_qdens[n]]<<std::endl;

    }
  }

  qdens.close();

}


void
Negf::parse_options(void)
{

  opt.pot_module = get_option("potential_simulation","none");

  ModelOptions& sol_opt = get_solver_options();

  opt.n_kT = sol_opt.get_option("n_kT", 10);

  opt.Emin = sol_opt.get_option("Emin",0.0);

  opt.Emax = sol_opt.get_option("Emax",0.0);

  opt.Estep = sol_opt.get_option("Estep",0.1);

  opt.Np_n.resize(2, 0);

  sol_opt.get_option("Np_contour", opt.Np_n);

  opt.delta = sol_opt.get_option("delta",2e-3);

  double deltaE =  sol_opt.get_option("deltaE",opt.delta/2.0);

  int Np = 2 * opt.n_kT/deltaE;

  opt.Np_real = sol_opt.get_option("Np_real",Np);

  opt.DEc = sol_opt.get_option("deltaEc", 0.25);

  opt.DEv = sol_opt.get_option("deltaEv", 0.25);

  opt.n_poles = sol_opt.get_option("Npoles", 0);

  opt.verbosity = sol_opt.get_option("verbosity",0);
  //sol_opt.check_unused();
  opt.writeLDOS = sol_opt.get_option("writeLDOS",false);

  opt.set_dirichlet_bc = get_option("dirichlet", false);

  if (!has_option("Np_real"))
  {
    std::cout<<"Np_real: "<<opt.Np_real<<std::endl;
  }
}

void
Negf::do_setup_solution_variables(void)
{


  declare_solution(ReorderPotential, REAL, NODES, "");
  declare_solution(elDensity, REAL, NODES, "1/cm^3");
  declare_solution(hDensity, REAL, NODES, "1/cm^3");
  declare_solution(eCurrentDensity, REAL, NODES, "A/cm^2");
  declare_solution(hCurrentDensity, REAL, NODES, "A/cm^2");

  if (plot_solution("Current"))
  {

    std::cout<<"Setup contact currents"<<std::endl;

    std::set<std::string> cnt_set;

    SimulationEnvironment::BoundaryIterator it = get_environment().boundaries_begin();
    const SimulationEnvironment::BoundaryIterator end = get_environment().boundaries_end();

    for ( ; it != end; ++it)
    {
      const Boundary* bd = (*it);
      if (bd != NULL)
      {
        std::cout<<bd->get_name()<<std::endl;

        std::string name(bd->get_name() + ".current");

        add_plot_variable(name);

        cnt_set.insert(name);
      }

      // now we declare them
      unsigned int dim = get_mesh().mesh_dimension();
      std::string units("A");
      if (dim == 1)
        units = "A/cm^2";
      else if (dim == 2)
        units = "A/cm";

      unsigned int id = static_cast<ID>(ContactCurrent);
      for (std::set<std::string>::iterator i(cnt_set.begin()); i != cnt_set.end(); ++i)
      {
        ++id;
        declare_solution_ext(*i, id, SolutionDescriptor::REAL,
            SolutionDescriptor::GLOBAL, units);
      }

    }

  }

}

void
Negf::plot_globaldata (void)
{
  if ( plot_solution("Current") )
  {
    std::string out_file = get_output_directory()+"/current.dat";

    std::fstream ff(out_file.c_str(), std::fstream::out);

    std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
    const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
    std::cout<<"contact"<<"       current:"<<std::endl;
    for (; it != end; ++it)
    {
      std::cout<<it->second->get_name()<<"  "<<_contact_current[it->second]<<std::endl;
      ff<< _contact_current[it->second] << "  ";
    }
    ff << std::endl;
    ff.close();
  }
}


void
Negf::ham_assemble(EquationSystems& es, const std::string& system_name)
{
  static_this->do_ham_assemble(es, system_name);
}


void
Negf::do_ham_assemble(EquationSystems& es, const std::string& system_name)
{
  SimulationInterface* model;
  ID sol_id;

  if (opt.pot_module != "none")
  {
    model = SimulationInterface::find_simulation(opt.pot_module);
    if (!model->is_solved() )
      throw SolveFailedException("Simulation "+opt.pot_module+" must be solved first");
  }
  else
  {
    model = NULL;
  }

  const double newconst = 0.5 * Hartree * bohr_radius/get_mesh_units() * bohr_radius/get_mesh_units();

  TensorValue<double> invMass(0.0);

  std::vector<double> V;

  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map = _sys_H->get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_face(FEBase::build(dim,fe_type));

  QGauss qrule(dim, THIRD);
  QGauss qrule_face(dim-1,FIRST);
  //QTrap qrule(dim);

  fe->attach_quadrature_rule(&qrule);
  fe_face->attach_quadrature_rule(&qrule_face);

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<Point>& face_normals = fe_face->get_normals();

  DenseMatrix<Number> Hr; // Interaction hamiltonian matrix real part
  DenseMatrix<Number> Sr; // Overlap matrix real part
  DenseMatrix<Number> Hi; // Interaction hamiltonian matrix immaginary part
  DenseMatrix<Number> Si; // Overlap matrix immaginary part

  std::vector<unsigned int> dof_indices,new_dof_indices;

  std::map<ID, QuantumContact*>::iterator qc_end = _quantum_contacts.end();

  //ACTIVATE QC
  activate_quantum_contacts();

  // Zero H and S
  _sys_H->matrix->zero();
  _sys_H->get_matrix("Hi").zero();
  _sys_S->matrix->zero();
  _sys_S->get_matrix("Si").zero();

  //ITERATION OVER ACTIVE DEVICE REGION
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;


    //-------------------------------------------------------
    //get effective mass tensor for elem
    NegfModel* negfmod = get_bulk_model<NegfModel>(elem);

    for (unsigned int band=0; band < negfmod->get_n_bands(); band++)
    {
      dof_map.dof_indices(elem, dof_indices, band);

      invMass = negfmod->get_inv_mass(band);

      // gets from model which band ("Ec" or "Ev")
      sol_id = model->get_solution_id(negfmod->get_band(band));

      const unsigned int n_dofs = dof_indices.size();

      fe->reinit(elem);

      //get potential from dd model
      V.resize(qrule.n_points());
      V.assign(qrule.n_points(), 0.0);

      std::map<ID, QuantumContact*>::iterator qc_it;

      if (model != NULL)
      {
        if ( (qc_it = _quantum_contacts.find(elem->subdomain_id())) != qc_end)
        {
          for (unsigned int qp=0; qp<q_point.size(); qp++)
          {
            std::pair<const Elem*, Point> pair = qc_it->second->project_on_boundary(elem, q_point[qp]);
            model->get_solution(pair.first, sol_id, V[qp], pair.second);
          }
        }
        else
        {
          model->get_solution(elem, sol_id, V, q_point);
        }
      }

      Hr.resize(n_dofs, n_dofs);
      Sr.resize(n_dofs, n_dofs);
      Hr.zero();
      Sr.zero();
      Hi.resize(n_dofs, n_dofs);
      Si.resize(n_dofs, n_dofs);
      Hi.zero();
      Si.zero();

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<phi.size(); i++)
          for (unsigned int j=0; j<phi.size(); j++)
          {
            Sr(i,j) += JxW[qp]* phi[i][qp] * phi[j][qp];
            Hr(i,j) += JxW[qp]* newconst * dphi[i][qp] * (invMass*dphi[j][qp]);
            Hr(i,j) += JxW[qp] * V[qp] * phi[i][qp] * phi[j][qp];
            Hr(i,j) += JxW[qp] * newconst * (_k_vec * (invMass * _k_vec)) * phi[i][qp] * phi[j][qp];
          }

      // DIRICHLET BC (Make sense only in 2D and 3D)---------------------------------
      if (dim>1)
      {
        unsigned int n_sides = elem->n_sides();


        for (short side = 0; side < n_sides; side++)
        {

          const ElementSide elside(elem->top_parent(), side);
          bool set_dirichlet = false;

          std::map<ID, QuantumContact*>::iterator it;

          if ( (it = _quantum_contacts.find(elem->subdomain_id())) != _quantum_contacts.end() )  //in quantum contact
          {
            if (elem->neighbor(side)==NULL)
            {
              set_dirichlet = opt.set_dirichlet_bc; //set to default value
              fe_face->reinit(elem, side);
              const Boundary* bd = get_boundary(it->second);
              std::string bc = bd->get_options().get_option("dirichlet","");
              if (bc!="")
              {
                int sign = (bc[0]=='-' ? -1: 1);
                int dir;
                if (bc[1]=='x') dir = 0;
                if (bc[1]=='y') dir = 1;
                if (bc[1]=='z') dir = 2;
                if  (sign*face_normals[0](dir) > 0) set_dirichlet = true;
              }
            }
            //std::cout<<"in qc "<<it->second->get_name()<<" "<<set_dirichlet<<std::endl;
          }

          /*if(_quantum_contacts.count(elem->subdomain_id())) //in quantum contact
          {
            if(elem->neighbor(side)==NULL)
            {
               //int sbd= (int) elem->subdomain_id();
               //std::cout<<"elem "<<sbd<<" : "<<elem->id()<<"  side:  "<<side<<std::endl;
               set_dirichlet = opt.set_dirichlet_bc; //set to default
            }
          }*/
          else if(_env->is_outer_boundary(elside))  //in device region
          {
            Boundary* bd = _env->get_boundary(elside);
            // bd==NULL => not defined in input (!quantum_contact)
            if ((bd == NULL && opt.set_dirichlet_bc) || _dirichlet_boundaries.count(bd) )
            {
              set_dirichlet = true;
            }
           }

          if (set_dirichlet)
          {
            //std::cout<<elem->id()<<"  ";
            //if (elem->neighbor(side) !=NULL) std::cout<<(elem->neighbor(side))->id()<<std::endl;
            //else std::cout<<std::endl;
            for (unsigned int nd = 0; nd < elem->n_nodes(); nd++)
            {
              if (elem->is_node_on_side(nd, side))
              {
                Hr(nd,nd) = newconst/0.01;
              }
            }
          }

        }

      }
      // -----------------------------------------------------------------------------------

      new_dof_indices.resize(n_dofs);


      for (unsigned int i=0; i< n_dofs; i++)
        new_dof_indices[i] = _inv_perm[dof_indices[i]];


      _sys_S->matrix->add_matrix(Sr, new_dof_indices);

      _sys_H->matrix->add_matrix(Hr, new_dof_indices);

      _sys_S->get_matrix("Si").add_matrix(Si, new_dof_indices);

      _sys_H->get_matrix("Hi").add_matrix(Hi, new_dof_indices);

    }//BAND LOOP
  }//ELEM LOOP

}

void
Negf::print_ham(std::string form)
{
  std::string outpath = SimulationOptions::scratch_path;

  if (form=="matlab")
  {
    _sys_H->matrix->print_matlab(outpath+"/Hr.m");
    _sys_S->matrix->print_matlab(outpath+"/Sr.m");
    _sys_H->get_matrix("Hi").print_matlab(outpath+"/Hi.m");
    _sys_S->get_matrix("Si").print_matlab(outpath+"/Si.m");
    //std::cout<<"print Matlab matrices"<<std::endl;
  }

}

void
Negf::print_Lib(void)
{
  ModelOptions& sol_opt = get_solver_options();
  double mu_n;
  double mu_p;
  double Ec = get_band_edge("Ec");
  double Ev = get_band_edge("Ev");
  double DeltaEc = -opt.DEc;
  double DeltaEv =  opt.DEv;
  unsigned int n_vars = _sys_H->n_vars();

  mu_n = mu_p = 0.0;     // set to 0.0 for libnegf:
                         // we set  Ef  as   Ef - mu (see below)

  double kbT = SimulationOptions::temperature * Constants::kb;
  double wght = 1.0;


  std::vector <double> Np_p(2);
  for (unsigned int i = 0; i < 2; i++)
  {
    Np_p[i] = 0.0;
  }

  unsigned int spin = 2;
  unsigned int nLDOS = 0;

  std::vector <unsigned int> LDOS(2*nLDOS);

  if (nLDOS > 0)
  {
    LDOS[0]=1; LDOS[1]=_device_n_dofs;
  }

  // phi: potential at boundaries (quantum contacts)
  // mu : electrochemical potential at boundaries (qc)
  std::vector <double> phi(_quantum_contacts.size(), 0.0);
  std::vector <double> mu(_quantum_contacts.size(), 0.0);


  mumin=10000;
  mumax=-10000;
  ID id = 0;
  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for (; it != end; ++it)
  {
    get_boundary_potentials(it->second, phi[id], mu[id]);

    if(mu[id]>mumax){mumax = mu[id];}
    if(mu[id]<mumin){mumin = mu[id];}
    _contact_potential[it->second] = mu[id];

    id++;
  }

  opt.Emin = sol_opt.get_option("Emin",-mumax-opt.n_kT*kbT);
  opt.Emax = sol_opt.get_option("Emax",-mumin+opt.n_kT*kbT);

  std::string outpath = SimulationOptions::scratch_path;
  std::string out_file = "negf.in";
  std::fstream ff(out_file.c_str(),std::fstream::out);

  ff<<"'"+outpath+"/Hr.m'"<<std::endl;
  ff<<"'"+outpath+"/Hi.m'"<<std::endl;
  ff<<"'"+outpath+"/Sr.m'"<<std::endl;
  ff<<"'"+outpath+"/Si.m'"<<std::endl;

  ff<<_quantum_contacts.size()<<std::endl;

  // Nblocks (if 0 are partitioned by the library)
  ff<<0<<std::endl;
  //ff<<floor(_device_n_dofs*n_vars/2)<<" "<<_device_n_dofs*n_vars<<std::endl;

  // contact end dofs:
  for (unsigned int i = 0; i <_quantum_contacts.size(); i++)
    ff<<_qc_n_dofs[i]*n_vars<<" ";
  ff<<std::endl;

  ff<<_device_n_dofs*n_vars<<" ";
  for (unsigned int i = 0; i <(_quantum_contacts.size()-1); i++)
    ff<<_qc_n_dofs[i]*n_vars<<" ";
  ff<<std::endl;

  ff<<mu_n<<" "<<mu_p<<std::endl;
  ff<<Ec+DeltaEc<<" "<<Ev+DeltaEv<<std::endl;
  ff<<0.0<<" "<<0.0<<std::endl;                // DEc DEv = 0 (are set above)
  ff<<opt.Emin<<" "<<opt.Emax<<" "<<opt.Estep<<std::endl;
  ff<<kbT<<std::endl;
  ff<<wght<<std::endl;

  for (unsigned int i = 0; i < 2; i++)
    ff<<opt.Np_n[i]<<" ";
  ff<<std::endl;

  for (unsigned int i = 0; i < 2; i++)
    ff<<Np_p[i]<<" ";
  ff<<std::endl;

  ff<<opt.Np_real<<std::endl;
  ff<<opt.n_kT<<std::endl;
  ff<<opt.n_poles<<std::endl;
  ff<<spin<<std::endl;
  ff<<opt.delta<<std::endl;
  ff<<nLDOS<<std::endl;

  for (unsigned int i = 0; i < 2*nLDOS; i++)
    ff<<LDOS[i]<<" ";

  ff<<std::endl;

  // Sets 0 for the Fermi level for LIBNEGF
  for (unsigned int i = 0; i < _quantum_contacts.size(); i++)
    ff<<0.0<<" ";
  ff<<std::endl;

  // takes -mu (ELECTROCHEMICAL potential) as LibNEGF will turn back the sign
  // Ef[k]-mu[k] = 0 - mu[k]   for LibNEGF
  for (unsigned int i = 0; i < _quantum_contacts.size(); i++)
    ff<<-mu[i]<<" ";
  ff<<std::endl;

  ff.close();

}

void
Negf::get_solution_secure(const Elem *elem, std::map<ID, std::vector<double>> &values,
    const std::vector<Point> &p)
{

  unsigned int np = p.size();
  const unsigned int dim = get_mesh().mesh_dimension();

  if (values.count(ReorderPotential))
  {
    const NumericVector<Number>& solution = _sys->get_solution_vector();
    const DofMap& dof_map = _sys->get_dof_map();
    ID u_var = _sys->variable_number("u0");
    FEType fe_type = _sys->variable_type(u_var);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &p);

    std::vector<unsigned int> dof_indices_u;

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    const unsigned int n_dofs = dof_indices_u.size();

    for (unsigned int n = 0; n < np; n++)
    {
      double u  = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        u += phi[i][n] * solution(dof_indices_u[i]); ;
      }

      values[ReorderPotential][n] = u;
    }
  }

  if (values.count(elDensity))
  {
    //std::cout<<"(Negf) Getting electron density"<<std::endl;
    _sys_H = &get_equation_system<TiberLinearSystem>(0);
    NumericVector<Number>& qdens = *_qdens_sys->solution;

    FEType fe_type = _qdens_sys->variable_type(0);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &p);

    DofMap& dof_map = _qdens_sys->get_dof_map();
    std::vector<unsigned int> dof_indices;
    dof_map.dof_indices(elem, dof_indices, 0);
    unsigned int n_dofs = phi.size();

    // the phi^2 factor comes from the fact that the more correct interpolation is the square
    // of the basis function, because the probability densities are the square of the states
    // NOTE: maybe one should check if this gives really a better result
    // NOTE: 2011-12-01 the above turned out to be wrong
    //  phi[i] rho[i][j] phi[j] ~~  phi[i] rho[i][i] phi[i]
    // rho[i][i] is the correct density at node i: phi is used linear interpolation.
    if (_k_int_density!=NULL)
    {

      for (unsigned int n = 0; n < np; n++)
      {
        double value = 0;

        for (unsigned int i = 0; i < n_dofs; i++)
        {
          value += phi[i][n] * density[dof_indices[i]];
        }

        values[elDensity][n] = value;
      }
    }
    else
    {
      for (unsigned int n = 0; n < np; n++)
      {
        double value = 0;

        for (unsigned int i = 0; i < n_dofs; i++)
          value += phi[i][n]  * qdens(dof_indices[i]);

        values[elDensity][n] = value;
      }
    }
  }

  if (values.count(eCurrentDensity))
  {
    //std::cout<<"(Negf) Getting e current density"<<std::endl;
    //given an element/point decide which contact it belongs to
    //values[eCurrentDensity] = std::vector<double>(np*3, 0.0);

    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      const ElementSide elside(elem, side);
      Boundary* b = _env->get_boundary(elside);

      if (b != NULL)
      {
        std::map<const Boundary*, QuantumContact*>::iterator it(_qc_boundaries.find(b));
        if (it !=  _qc_boundaries.end())
        {
          unsigned int dummy;
          Point point_current = _contact_current[it->second] * it->second->get_normal(dummy);

          // assign the same current to all points
          for (unsigned int n = 0; n < np; n++)
          {
            values[eCurrentDensity][n*3]   =  point_current(0);
            values[eCurrentDensity][n*3+1] =  point_current(1);
            values[eCurrentDensity][n*3+2] =  point_current(2);
          }
        }
        break;
      }

    }

  }

}

void
Negf::get_solution_secure(std::map<ID, std::vector<double> >& values)
{

  std::vector<std::string> tokens;

  std::map<ID, std::vector<double> >::iterator mapit(values.begin());
  const std::map<ID, std::vector<double> >::iterator mapend(values.end());

  for ( ; mapit != mapend; ++mapit)
  {

    ID id = mapit->first;

    const SolutionDescriptor& descr = get_solution_descriptor(id);
    Utils::tokenize(descr.name(), tokens);

    std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
    const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
    for (; it != end; ++it)
    {
      if (tokens[0]==it->second->get_name())
      {
        values[id] = std::vector<double>(1,  _contact_current[it->second]);
        break;
      }
    }
  }

}

// ACTIVATE TEMPORARILY ALL ELEMENTS IN QUANTUM CONTACTS
void
Negf::activate_quantum_contacts(void)
{
  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for (;it != end;++it)
    it->second->activate_elements();
}
// DEACTIVATE TEMPORARILY ALL ELEMENTS IN QUANTUM CONTACTS
void
Negf::deactivate_quantum_contacts(void)
{
  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for (;it != end;++it)
    it->second->inactivate_elements();
}
// -----------------------------------------------------------------------
// Assemble a Laplace problem in active regions and use solution to
// reorder dofs
void
Negf::reorder(void)
{
  //std::cout<<"Reorder dofs"<<std::endl;

  MeshBase& mesh = get_mesh();

  unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_var_reord = 1;
  //std::<<"n vars "<<n_vars<<std::endl;

  std::vector<ID> u_vars(n_var_reord,0);

  for (unsigned int n=0; n< n_var_reord; n++)
  {
    std::ostringstream var_str;
    var_str << "u" << n;
    std::string name = var_str.str();
    _sys->add_variable(name, FIRST, LAGRANGE);
    u_vars[n] = _sys->variable_number(name);
  }

  _sys->attach_assemble_function(reorder_assemble);

  _sys->init();

  _sys->solve();

  //std::cout<<"Laplace solved"<<std::endl;

  const NumericVector<Number>& solution = _sys->get_solution_vector();

  const DofMap& dof_map = _sys->get_dof_map();

  std::vector<unsigned int> dof_indices_u;

  unsigned int sol_size = solution.size();

  // setup initial permutation vector as identitiy
  // the vector runs only on the device region where the
  // dof reordering is performed
  _perm.resize(_device_n_dofs, 0);
  for (unsigned int i = 0; i < _device_n_dofs; i++)
    _perm[i]=i;

  // Number of dofs of the first band
  unsigned int n_dofs = _qc_n_dofs[_quantum_contacts.size()-1];

  //std::cout<<"Sorting"<<std::endl;
  //std::cerr<<"n_dofs: "<<n_dofs<<std::endl;
  //std::cerr<<"dev: "<<_device_n_dofs<<std::endl;

  std::sort(_perm.begin(), _perm.end(), compare);

  // ========================================================================
  // Initial Dofs:
  //
  // ---------device----------||-cont--||--------device--------||---cont---
  // |  1 |  2 |  3 |  4 |  5 || 6 | 7 || 8 | 9 | 10 | 11 | 12 || 13 | 14 |
  //
  // After Sorting (which is done only on device):
  //
  // | 4  |  1 |  2 |  5 |  3 |
  // (P[1]=4; P[2]=1; P[3]=2; P[4]=5; P[5]=3)
  //
  // We need to invert the permutation _P[P[i]] = i
  // _P[4] = 1; _P[1] = 2;  _P[2] = 3; _P[5] = 4; _P[3] = 5;
  //
  // |_P[4]| _P[1]| _P[2]| _P[5]| _P[3]|
  //
  // FOR MORE VARIABLES AND DOFs
  //
  // ---------device----------||-cont--||--------device--------||---cont---
  // |  1 |  2 |  3 |  4 |  5 || 6 | 7 || 8 | 9 | 10 | 11 | 12 || 13 | 14 |
  //
  // | 4  | 11 |  1 |  8 |  2 | 9  | 5  | 12 | 3 | 10 || 6 | 13 | 7  | 14 |
  //
  // =========================================================================

  // Invert permutation to compute _inv_perm
  _inv_perm.clear();
  _inv_perm.resize(n_dofs*n_vars, 0);

  for (unsigned int i = 0; i < _device_n_dofs ; i++)
    for (unsigned int k = 0; k < n_vars; k++)
      _inv_perm[_perm[i]+k*n_dofs]= i*n_vars+k;

  // Reset dofs in QC as an identity
  for (unsigned int i = 0; i < n_dofs - _device_n_dofs ; i++)
    for (unsigned int k = 0; k < n_vars; k++)
      _inv_perm[_device_n_dofs+i+k*n_dofs] = _device_n_dofs*n_vars + i*n_vars + k;

  // Print permutation data
  activate_quantum_contacts();
  ID assign = 0;
  std::vector<unsigned int> temp(sol_size,sol_size);

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices_u);

    const unsigned int n_dofs = dof_indices_u.size();

    for (unsigned int i = 0; i < n_dofs; i++)
    {
      if (temp[dof_indices_u[i]] == sol_size)
      {
        temp[dof_indices_u[i]] = assign;
        assign++;
      }
    }
  }
}

void
Negf::reorder_assemble(EquationSystems& es, const std::string& system_name)
{
  static_this->do_reorder_assemble(es, system_name);
}

//Approch to Laplace problem
void
Negf::do_reorder_assemble(EquationSystems& es, const std::string& system_name)
{
  _device_n_dofs = 0;
  _qc_n_dofs.resize(_quantum_contacts.size(), 0);

  std::map<const Boundary*, int> boundary_ids;
  int id = 0;
  std::map<const Boundary*, QuantumContact*>::iterator it = _qc_boundaries.begin();
  const std::map<const Boundary*, QuantumContact*>::iterator end = _qc_boundaries.end();
  for( ; it != end; ++it, ++id)
  {
    boundary_ids[it->first] = id;
  }


  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map = _sys->get_dof_map();

  unsigned int n_vars = _sys->n_vars();

  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

  QGauss qrule(dim, THIRD);

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  const double penalty = 1.e10;

  // DEACTIVATE TEMPORARILY ALL ELEMENTS IN QUANTUM CONTACTS
  deactivate_quantum_contacts();
  for (unsigned int k=0; k < n_vars; k++)
  {
    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    _dev.resize(1);

    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      fe->reinit(elem);

      dof_map.dof_indices(elem, dof_indices, k);
      const unsigned int n_dofs = dof_indices.size();

      Ke.resize(n_dofs, n_dofs);
      Fe.resize(n_dofs);
      Ke.zero();

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<phi.size(); i++)
          for (unsigned int j=0; j<phi.size(); j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      Fe.zero();

      // Set Boundary conditions at contact boundaries
      // Currently only 2 contacts makes sense
      for (unsigned int side=0; side<elem->n_sides(); side++)
      {
        const ElementSide elside(elem->top_parent(), side);

        Boundary* bd = _env->get_boundary(elside);

        if (bd != NULL && _qc_boundaries.count(bd))
        {
          ID cc = boundary_ids[bd];

          for(unsigned int n = 0; n< n_dofs; ++n)
          {
            if (elem->is_node_on_side(n,side))
            {
              for (unsigned int nc = 0; nc < n_dofs; nc++)
                Ke(n,nc) = 0.0;

              Ke(n,n) = penalty;

              Fe(n) = penalty*cc; //fix potential
            }
          }
        }
      }

      _sys->matrix->add_matrix(Ke, dof_indices);
      _sys->rhs->add_vector(Fe, dof_indices);

      // Update number of dofs in device region
      for (unsigned int n = 0; n <n_dofs; n++)
        if (_device_n_dofs < dof_indices[n]+1) _device_n_dofs = dof_indices[n]+1;

    }
  }

  // Iterates over INACTIVE QuantumContact regions
  {
    for (unsigned int k=0; k < n_vars; k++)
    {
      MeshBase::const_element_iterator       el     = mesh.not_active_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.not_active_elements_end();

      _qc.resize(1);
      for ( ; el != end_el; ++el)
      {
        const Elem* elem = *el;

        ID sub_id = elem->subdomain_id();
        // Checks whether the element is within a QuantumContact
        if (_quantum_contacts.count(sub_id))
        {

          dof_map.dof_indices(elem, dof_indices, k);
          const unsigned int n_dofs= dof_indices.size();

          fe->reinit(elem);
          Ke.resize(n_dofs, n_dofs);
          Fe.resize(n_dofs);
          Ke.zero();
          Fe.zero();

          std::string qc_name = _device->get_region_name(sub_id);

          unsigned int cc = boundary_ids[_env->get_boundary(qc_name)];  // cc = 0 or 1

          for(unsigned int n = 0; n< n_dofs; ++n)
          {
            Ke(n,n) = penalty;
            Fe(n) = penalty * cc;

            // Update number of dofs in quantum contact. It begin from number of dofs in device region
            if (_qc_n_dofs[cc] <= dof_indices[n]+1) _qc_n_dofs[cc] = dof_indices[n]+1;
          }

          _sys->matrix->add_matrix(Ke, dof_indices);
          _sys->rhs->add_vector(Fe, dof_indices);

        }

      }

      // Reorder the qc_n_dofs such that they are with increasing order
      for (unsigned int i = 0; i <(_quantum_contacts.size()-1); i++)
        if (_qc_n_dofs[i]>_qc_n_dofs[i+1])
          std::swap(_qc_n_dofs[i],_qc_n_dofs[i+1]);

    }
  }
}

//Need to compare two solution to reorder dof indices
bool
Negf::compare(ID i, ID j)
{
  return static_this->do_compare(i, j);
}

bool
Negf::do_compare(ID i, ID j)
{
  const NumericVector<Number>& solution = _sys->get_solution_vector();
  return (solution(i) < solution(j));
}

void
Negf::get_boundary_potentials(QuantumContact* qc, double& av_V, double& av_mu)
{
  PotentialInterface model;
  av_V = 0.0;
  av_mu = 0.0;

  if (opt.pot_module == "none") return;
  else model.set_simulation(opt.pot_module);

  std::vector<double> V;

  MeshBase& mesh = get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  AutoPtr<FEBase> fe( FEBase::build(dim, FEType() ));

  QGauss qrule(dim, FIRST); // Order 0 rule because in this way we take centroid's normal

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Point>& q_point = fe->get_xyz();

  ID qc_id =qc->get_id();

  double volume = 0.0;
  MeshBase::const_element_iterator it= mesh.active_elements_begin();
  MeshBase::const_element_iterator it_end  = mesh.active_elements_end();
  for ( ; it != it_end ; ++it) //loop over k space elements
  {
    const Elem* elem = *it;
    if (elem->subdomain_id() == qc_id)
    {
      fe->reinit(elem);
      for (unsigned int qp=0; qp<q_point.size(); qp++)
      {
        volume +=  elem->volume();
      }
    }
  }

  it = mesh.active_elements_begin();
  it_end  = mesh.active_elements_end();

  for ( ; it != it_end ; ++it) //loop over k space elements
  {
    const Elem* elem = *it;

    if (elem->subdomain_id() == qc_id)
    {
      fe->reinit(elem);

      for (unsigned int qp=0; qp<q_point.size(); qp++)
      {
        std::pair<const Elem*, Point> pair = qc->project_on_boundary(elem, q_point[qp]);

        av_V += model.get_potential(pair.first, pair.second) * elem->volume()/volume ;
        av_mu += model.get_el_chem_potential(pair.first, pair.second) * elem->volume()/volume;
      }
    }
  }
}


double Negf::get_band_edge(const std::string& band) const
{
  SimulationInterface* model;
  MeshBase& mesh = get_mesh();

  if (opt.pot_module != "none")
  {
    model = SimulationInterface::find_simulation(opt.pot_module);
    if (!model->is_solved() )
      throw SolveFailedException("Simulation "+opt.pot_module+" must be solved first");
  }
  else
  {
    model = NULL;
  }

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  assert(el != end_el);

  double band_edge = get_band_edge(model, band, *el);
  ++el;

  for (; el != end_el ; ++el )
  {
    const Elem* elem = *el;

    double temp = get_band_edge(model, band, elem);

    if (band=="Ec")
    {
      if (band_edge > temp)
        band_edge = temp;
    }
    else
    {
      if (band_edge < temp)
        band_edge = temp;
    }
  }

  return (band_edge);

}

inline double Negf::get_band_edge(SimulationInterface* model, const std::string& band, const Elem* elem) const
{
  std::vector<double> values(elem->n_nodes());
  std::vector<Point> p(elem->n_nodes());

  for (unsigned int i = 0; i < elem->n_nodes(); ++i)
    p[i] = elem->point(i);

  ID band_edge_ID = model->get_solution_id(band);

  model->get_solution(elem, band_edge_ID, values, p);

  double bedge = values[0];

  for (size_t i = 1; i < elem->n_nodes(); ++i)
  {
    double temp = values[i];
    if(band == "Ec")
      bedge = (temp < bedge) ? temp : bedge;
    else
      bedge = (temp > bedge) ? temp : bedge;
  }

  return bedge;
}





void
Negf::apply_dirichlet_bc(void)
{
  const MeshBase* mesh = &get_mesh();

  unsigned int dim = mesh->mesh_dimension();

  MeshBase::const_element_iterator it = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh->active_local_elements_end();

  //activate_contacts(); // dirichelet must go on contacts as well.

  dirichlet_dofs.clear();

  DofMap& dof_map = _sys_H->get_dof_map();

  unsigned int n_var = dof_map.n_variables();

  std::vector<unsigned int> dof_indices;


  // iterates over all elements and mark the external sides
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    unsigned int n_sides;

    if ( dim > 1 )
      n_sides = elem->n_sides();
    else
      n_sides = elem->n_nodes();


    for (short side = 0; side < n_sides; side++)
    {

      Elem* el1 = elem->neighbor(side);
      bool side_is_external = false;

      // mesh refinements are not handled
      if (el1 == NULL)
      {
         side_is_external = true;
      }
      else  if ( !(el1->active()) )
      {
         side_is_external = true;
      }

      if (side_is_external)
      {

        if (dim > 1)
        {//2D/3D
          for (unsigned int nd = 0; nd < elem->n_nodes(); nd++)
          {
            if (elem->is_node_on_side(nd, side))
            {
              const Point& node = elem->point(nd);
              for (short band = 0 ; band <  n_var; band++)
              {
                dof_map.dof_indices(elem, dof_indices, band);
                dirichlet_dofs.insert(dof_indices[nd]);
              }
            }
          }

        }
        else
        {//1D
          for (short band = 0 ; band < n_var; band++)
          {
            dof_map.dof_indices (elem, dof_indices, band);
            dirichlet_dofs.insert(dof_indices[side]);
          }
        }

      }

    }// close side loop

  }// close element loop
}


const Boundary*
Negf::get_boundary(const QuantumContact* qc)
{
  std::map<const QuantumContact*, const Boundary*>::iterator it = _bd_map.find(qc);
  return it->second;
}
