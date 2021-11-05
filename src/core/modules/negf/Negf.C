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

#include "libnegf/NegfWrapper.h"

// To be able to compile as module
#include "TiberModule.h"

// Basic include files needed for the mesh functionality.
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
// Define generic quadrature rules.
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/mesh.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/tensor_value.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"


// C++ includes
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cassert>
#include <cmath>

using namespace Constants;
using namespace libMesh;
using namespace std;

namespace {
  struct sortclass{
      sortclass(const std::vector<Atom>& atoms) : _atoms(atoms) {}
      ~sortclass(){};
      bool operator() (int i, int j) { return (_atoms[i].get_position()(0)<_atoms[j].get_position()(0)); }
      const std::vector<Atom>& _atoms;
  };
}


 


namespace
{
  double fermi(double Energy, double Fermi_energy, double Temperature)
  {
    double T_EV = Temperature * Constants::k_Boltzmann;
    double exp_arg =  (Energy - Fermi_energy)/T_EV;

    double occupation;

    if (exp_arg > 35)
      occupation = std::exp(-exp_arg);
    else
      occupation = 1.0/(std::exp(exp_arg) + 1.0);

    return occupation;
  }
}



Negf::Negf(const ModelOptions& options) :
  SimulationInterface(options),
  _reorder_assembly(this),
  _device_n_dofs(0),
  _k_int_density(NULL),
  _k_int_current(NULL),
  _ext_module(NULL)
{
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
  return(new Negf(options));
}

PhysicalModel*
Negf::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  if (!options.has_submodel("hamiltonian"))
  {
    Messages::error("Physics must contain hamiltonian model");  
  }
  return  NegfModel::create(mat, options);
}


//! The initialization
void
Negf::do_init(void)
{

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

        qc->set_neighbor_map();

      }
    }
  }

  // Now we assign an internal, successive integer to every contact. This
  // number will be used as boundary condition for the reordering. We do
  // this based on an ordering in coordinates.

  std::map<ID, Point> min_coord;
  std::vector<ID> qids;

  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for( ; it != end; ++it)
  {
    QuantumContact& qc = *it->second;
    ID qid = it->first;

    qids.push_back(qid);

    min_coord[qid] = Point(1e16, 1e16, 1e16);

    for (auto s = qc.contact_elements_begin(); s != qc.contact_elements_end(); ++s)
    {
      const Elem* el = s->first;
      for (int n = 0; n < el->n_nodes(); ++n)
      {
        if (el->point(n) < min_coord[qid])
          min_coord[qid] = el->point(n);
      }
    }
  }
  
  for(unsigned int i = 0; i < qids.size() - 1; ++i)
  {
    for(unsigned int j = i + 1; j < qids.size(); ++j)
    {
      if (min_coord[qids[j]] < min_coord[qids[i]])
        std::swap(qids[i], qids[j]);
    }
  }

  Messages m;
  m.info("Quantum contacts with IDs for reordering:");
  m.indent();
  unsigned int id = 0;
  for(unsigned int i = 0; i < qids.size(); ++i, ++id)
  {
    _bd_num[_bd_map[_quantum_contacts[qids[i]]]] = id;

    ostringstream os;
    os << _quantum_contacts[qids[i]]->get_name() << " -> " << i;
    m.info(os.str());
  }

  _qc_n_dofs.resize(_quantum_contacts.size(), 0);
  
  init_hamil();
  
}

void
Negf::init_hamil(void)
{

  // setup external module as Hamiltonian generator
  
  const MeshBase& mesh = get_mesh();
  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  NegfModel* negfmod;

  negfmod = get_bulk_model<NegfModel>(*el);

  _hamil_type = negfmod->get_model_name(0);

  std::cout<<"(negf) init: "<< _hamil_type  <<std::endl;

  std::string sim = negfmod->get_simulation(0);
  std::cout<<"(negf) sim: "<< sim <<std::endl;

  _ext_module = dynamic_cast<EigenvalueProblem*>(find_simulation(sim));

  if ( _hamil_type == "etb")  
  {
     init_etb_hamil();
  }
  else if ( _hamil_type == "efa")
  {
     init_efa_hamil();
  }  
  else 
  {
    Messages::error("undefined hamiltonian type "+_hamil_type);
  }
    
  
  std::cout<<"(negf) init done. " <<std::endl;
  
}


//! The initialization
void
Negf::init_efa_hamil(void)
{

  unsigned int n_bands;

  if (_ext_module == NULL)
  {
    throw InitFailedException("NEGF module needs an external"
        " provider of the Hamiltonian.");
  }
  else
  {
    if (!_ext_module->is_initialized())
      _ext_module->init();
    n_bands = _ext_module->get_number_of_bands();
  }

  // Setup a simple effective mass Hamiltonian
  std::cout<<"(negf) create eq_sys with "<<n_bands<<" bands"<<std::endl;

  ID id = create_equation_system("linear","");
  _sys_H = &get_equation_system<TiberLinearSystem>(id);

  // attach a variable for each subband
  for (unsigned int band=0; band < n_bands; band++)
  {
    std::stringstream out;
    out<<"phi"<<band;
    _sys_H->add_variable(out.str(), FIRST, LAGRANGE, &get_region_ids());
  }


  // attach system for reorder dofs
  std::cout<<"(negf) get eq sys for reorder"<<std::endl;
  id = create_equation_system("linear","");
  _sys = &get_equation_system<TiberLinearSystem>(id);

  // We add a second system just to contain the density
  if (plot_solution(elDensity) || plot_solution(hlDensity) ||
      plot_solution("LDOS"))
  {
    std::cout<<"(negf) create eq sys for elDensity"<<std::endl;
    id = create_equation_system("linear","");
    _qdens_sys = &get_equation_system<TiberLinearSystem>(id);
    _qdens_sys->add_variable("edens", libMeshEnums::FIRST, LAGRANGE, &get_region_ids());
    _qdens_sys->add_variable("hdens", libMeshEnums::FIRST, LAGRANGE, &get_region_ids());
    _qdens_sys->init();
  }

  //std::cout<<"(negf) init H and S"<<std::endl;

  _sys_H->init();

  std::cout<<"(negf) init k-integration"<<std::endl;
  init_k_space_integration();

  std::cout<<"(negf) reorder dofs"<<std::endl;
  reorder(); // dof indices reorder

}


void
Negf::init_etb_hamil(void)
{
   AtomisticStructure* as = _ext_module->get_atomistic_structure();
   unsigned int N_atoms = as->get_N_atoms();
   std::vector<Atom>& atoms = as->get_structure_atoms();

   //Reordering for negf. Currently works for 2 contacts
   
   unsigned int last_dev = -1;
   std::map<ID,unsigned int> last_cont;
   std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
   const std::map<ID, QuantumContact*>::iterator qc_end = _quantum_contacts.end();
   for (it; it!=qc_end; it++)
   {  
     last_cont.insert(std::pair<ID,unsigned int>(it->first,-1) );
   }

   _inv_perm.resize(N_atoms);

   // Count all atoms in device and contacts regions 
   // Set inv_perm in device region (order left unchanged)
   for (unsigned int i=0; i< N_atoms; i++)
   {
       it = _quantum_contacts.find(atoms[i].get_region_ID());         
       if (it == qc_end)
       {
          _inv_perm[i] = ++last_dev;
       }
       if (it != qc_end)
       {
          (last_cont[it->first])++;
       } 
   } 

   // Set inv_perm in contacts 
   // contact with normal>0 is left unchanged the other is reversed
   unsigned int count1 = 0;
   unsigned int count2 = N_atoms;
   for (unsigned int i=0; i< N_atoms; i++)
   {
       it = _quantum_contacts.find(atoms[i].get_region_ID());         
       if (it != qc_end)
       {
          if (it->second->get_normal()(0) > 0) 
          {
             _inv_perm[i]= last_dev+ (++count1); 
          }     
          else
          {
             _inv_perm[i]= (--count2);
          }
       }
   }

   // invert to find perm 
   _perm.resize(N_atoms);
   for (unsigned int i=0; i< N_atoms; i++)
   {
     _perm[_inv_perm[i]]=i;
   }

   sortclass sortobj(atoms);
   std::sort(_perm.begin(),_perm.begin()+last_dev+1,sortobj);

   //for (unsigned int i=0; i< N_atoms; i++)
   //{
   //   std::cout<<"atom "<<i<<" reg "<<atoms[i].get_region_ID()
   //            <<"iperm "<<_inv_perm[i]<<" perm "<<_perm[i]<<std::endl;
   //}
   
  
   // ------------------------------------------------------------------------
   // count assign dofs in all regions
   _device_n_dofs = last_dev+1;

   std::cout<<"Number of atoms: "<<N_atoms<<std::endl;
   std::cout<<"Device: "<<_device_n_dofs<<std::endl;
   
   it = _quantum_contacts.begin();
   for (it; it!=qc_end; it++)
   {  
     std::cout<<"Contact: "<<it->first<<"  "<<last_cont[it->first]+1<<std::endl;
   }
  
   _qc_n_dofs[0] =  _device_n_dofs;
   _qc_n_dofs[1] =  _device_n_dofs;

   it = _quantum_contacts.begin();
   for (it; it!=qc_end; it++)
   {  
     //std::cout<<"Normal: "<<(it->second)->get_normal()(0)<<std::endl;
     if (it->second->get_normal()(0) > 0)
     {
       _qc_n_dofs[0] = _qc_n_dofs[1] + last_cont[it->first]+1;
     }
     else
     {
       _qc_n_dofs[1] = _qc_n_dofs[0] + last_cont[it->first]+1;
     }
   }
   
   //std::cout<<"write: "<<std::endl;
   for (ID i=0; i<2; i++)
   {
      std::cout<<"Conts: "<<i+1<<"  "<<_qc_n_dofs[i]+1<<std::endl;
   }

  std::cout<<"(negf) init k-integration"<<std::endl;
  init_k_space_integration();
}


void
Negf::do_reinit(void)
{
  //std::cout<<"(negf) clean up libnegf"<<std::endl;
  //_libnegf->clean_libnegf();
}

void
Negf::init_k_space(ModelOptions& kopts)
{
    unsigned int dim = kopts.get_option("k_space_dimension", 0);
    unsigned int mesh_dim = get_mesh().mesh_dimension();

    if ( _hamil_type == "etb" ) 
    {
      const AtomisticStructure* as = _ext_module->get_atomistic_structure();
      // Define k-space from cartesian lattice vectors
      double cf = 0.1;  //conversion factor from Angstrom to nm
      const std::vector<libMesh::RealVectorValue>& vectors = as->get_lattice_vectors();
      std::vector<double> r1(3,0.0);
      std::vector<double> r2(3,0.0);
      switch (dim)
      {
         case 2:
           r1[0] = cf*vectors[1](0); r1[1] = cf*vectors[1](1); r1[2] = cf*vectors[1](2);
           r2[0] = cf*vectors[2](0); r2[1] = cf*vectors[2](1); r2[2] = cf*vectors[2](2);
           std::cout<<"(Negf) lattice vec:"<<r1[0]<<" "<<r1[1]<<" "<<r1[2]<<std::endl;
           std::cout<<"(Negf) lattice vec:"<<r2[0]<<" "<<r2[1]<<" "<<r2[2]<<std::endl;
           kopts.set_option("r1",r1);  
           kopts.set_option("r2",r2);  
           break;

         case 1:
           r1[0] = cf*vectors[2](0); r2[1] = cf*vectors[2](1); r2[2] = cf*vectors[2](2);
           kopts.set_option("r1",r1);  
           break;

         default:
           break;
      }
    }
    else
    {
      // it is efa, we setup just some lattice vectors

      // Modifiche Alex rimaste appese. Ricordo problema se mesh in um !
      //double cf = 1e-10/get_mesh_units(); //conversion from Angstrom to M.U.
      // the factor of 2 is due to the fact that the k-space limits
      // in Kspace are 1/2
      //cf *= 20;
      // these are the real space lattice vectors, 1 cf = 1 nm
      //RealVectorValue a(cf, 0, 0), b(0, cf, 0), c(0, 0, cf);

      // we use 1 nm^-1 as default k-space extension (max k)
      // these are the real space lattice vectors, in nm
      // why pi? Because then the default max k becomes 1 ( = 2*pi/(2*a) ), and
      // k_max can be interpreted in nm^-1
      RealVectorValue a(M_PI, 0.0, 0), b(0.0, M_PI, 0.0), c(0.0, 0.0, M_PI);

      auto bbox = get_environment().get_bounding_box();
      if (get_option("x-periodicity", false) && (mesh_dim > 0))
        a(0) = (bbox.second(0) - bbox.first(0)) * get_mesh_units() * 1e9;
      if (get_option("y-periodicity", false) && (mesh_dim > 1))
        b(1) = (bbox.second(1) - bbox.first(1)) * get_mesh_units() * 1e9;
      if (get_option("z-periodicity", false) && (mesh_dim > 2))
        c(2) = (bbox.second(2) - bbox.first(2)) * get_mesh_units() * 1e9;

      switch (dim)
      {
        case 1:
          if (mesh_dim == 3)
          {
            if (get_option("x-periodicity", false))
              c = a;
            else if (get_option("y-periodicity", false))
              c = b;
          }
          kopts.set_option("r1", c);
          break;

        case 2:
          if (mesh_dim == 3)
          {
            if (!get_option("z-periodicity", false))
            {
              c = b;
              b = a;
            }
            else if (get_option("x-periodicity", false))
            {
              b = a;
            }
          }
          else if (mesh_dim == 2)
          {
            if (get_option("x-periodicity", false))
              b = a;
          }
          kopts.set_option("r1", b);
          kopts.set_option("r2", c);
          break;

        case 3:
          kopts.set_option("r1", a);
          kopts.set_option("r2", b);
          kopts.set_option("r3", c);
          break;

        default:
          break;

      }
    }

}

void
Negf::init_k_space_integration(void)
{
  int dim = get_mesh().mesh_dimension();

  unsigned int k_dim = 3 - dim;
  if (get_option("x-periodicity", false))
    k_dim++;
  if (get_option("y-periodicity", false))
    k_dim++;
  if (get_option("z-periodicity", false))
    k_dim++;

  k_dim = min(k_dim, 3u);


  if (get_options().has_submodel("k_integration_density"))
  {
    //-----------------DENSITY ---------------------------------------------
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_density"));
    ModelOptions kopts;
    kopts = it->second;

    //if (it != get_options().submodels_end("k_integration_density"))

    //kopts.set_option("mesh_units", get_mesh_units());

    kopts.set_option("k_space_dimension",k_dim);
    if ( !kopts.find_option("verbose"))
    {
      kopts.set_option("verbose", SimulationOptions::verbose() );
    }

    init_k_space(kopts); 

    _k_int_density = KspaceIntegration::create(this, &Negf::calculate_for_k_point, kopts,
                                               get_communicator(),
                                               get_solver_communicator());

    if (_k_int_density == NULL)
      throw InitFailedException("Could not create k-integration");

    _k_int_density->init();
  }
  //------------------CURRENT --------------------------------------------------
  if (get_options().has_submodel("k_integration_current"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_current"));
    ModelOptions& kopts = it->second;

    //kopts.set_option("mesh_units", get_mesh_units());
    kopts.set_option("k_space_dimension", k_dim);
    if (! kopts.find_option("verbose"))
    {
      kopts.set_option("verbose", SimulationOptions::verbose() );
    }

    init_k_space(kopts); 
    
    _k_int_current = KspaceIntegration::create(this, &Negf::calculate_for_k_point, kopts,
                                               get_communicator(),
                                               get_solver_communicator());

    if (_k_int_current == NULL)
      throw InitFailedException("Could not create k-integration");

    _k_int_current->init();

  }
}

void
Negf::setup_hamil(void)
{
  if ( _hamil_type == "efa")
  {
      setup_efa_hamil();
  }
  else if ( _hamil_type == "etb")
  {
      setup_etb_hamil();
  }

  int nrow = _ext_module->get_H_dim();

  {
    int nnz = _ext_module->get_H_nnz();

    std::vector<int> IA(nrow+1,0);
    std::vector<int> JA(nnz,0);
    std::vector<Complex> A(nnz);

    Messages::info("Passing Hamiltonian ... ", 0);

    if ( _hamil_type == "etb" )
      _ext_module->get_H_csr(A, JA, IA);
    else
      _ext_module->get_H_csr(A, JA, IA, _perm);

    _libnegf->set_H_csr(nrow, A, JA, IA);

    Messages::info("done.");
  }

  if (_ext_module->is_generalized())
  {

    // here we do not have a method to obtain nnz exactly,
    // but it will be adjusted inside get_S_csr() anyway
    std::vector<Complex> S(nrow);
    std::vector<int> IS(nrow+1,0);
    std::vector<int> JS(nrow,0);

    Messages::info("Passing overlap ... ", 0);

    if ( _hamil_type == "etb" )
      _ext_module->get_S_csr(S, JS, IS);
    else
      _ext_module->get_S_csr(S, JS, IS, _perm);


    _libnegf->set_S_csr(nrow, S, JS, IS);

    Messages::info("done.");
  }
  else
  {
    _libnegf->set_S_id(nrow);
  }

  // we have to reinitialize some structures in libnegf
  if (get_option("print_matrices",false))
    _libnegf->print_mat();

  _libnegf->init_structure(_quantum_contacts.size(),
                           _surfstart, _surfend, _contend,
                           _end_blocks.size(), _end_blocks);
}


void
Negf::setup_efa_hamil(void)
{
  do_reinit();
   
  Point k_point; 
  for(short i=0;i<3;i++) k_point(i) = _k_vec(i);

  _ext_module->set_k_point(k_point);

  _ext_module->reinit();
 
  _ext_module->assemble();
 
}


void
Negf::setup_etb_hamil(void)
{
  do_reinit();
   
  //apply permutation to atomistic structure
   _ext_module->get_atomistic_structure()->reorder(_perm);
  
  Point k_point; 
  for(short i=0;i<3;i++) k_point(i) = _k_vec(i);

  _ext_module->set_k_point(k_point);

  _ext_module->reinit();

  _ext_module->assemble();

  //_ext_module->print_H(get_scratch_directory());

}



void
Negf::setup_negf(void)
{

  // first time it will init
  _libnegf->init();

  _libnegf->set_scratch_path(get_scratch_directory());

  _libnegf->set_output_path(get_output_directory());

  // initialize the contacts
  _libnegf->init_contacts(_quantum_contacts.size());

  // new way via memory
  // first get the parameter strcture with defaults
  NegfWrapper::Parameters params;
  _libnegf->get_parameters(params);

  // now set values
  params.verbose = opt.verbosity;
  //params.readOldDM_SGFs = ;
  //params.readOldT_SGFs = ;
  params.g_spin = _ext_module->get_degeneracy();
  params.delta = opt.delta;
  //params.deltaModel = ;
  //params.wmax = ;
  //params.dos_delta = ;

  //std::vector <double> phi(_quantum_contacts.size(), 0.0);
  //std::vector <double> mu_n(_quantum_contacts.size(), 0.0);
  //std::vector <double> mu_p(_quantum_contacts.size(), 0.0);

  double kbT = SimulationOptions::temperature * Constants::kb;
  ModelOptions& sol_opt = get_solver_options();

  mumin=10000;
  mumax=-10000;
  ID id = 0;
  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for (; it != end; ++it)
  {
    
    if (id >= MAXNCONT)
    {
      ostringstream os;
      os << "Maximum number of contacts in NEGF is limited to "
        << MAXNCONT;
      throw InitFailedException(os.str());
    }

    double phi, mu_n, mu_p;
    get_boundary_potentials(it->second, phi, mu_n, mu_p);

    params.mu_n[id] =  mu_n;
    params.mu_p[id] = -mu_p;
    //params.mu[id] = ; // not needed, for DFTB

    // for now T is the same everywhere
    params.kbt_dm[id] = kbT;
    params.kbt_t[id]  = kbT;

    if (mu_n > mumax) mumax = mu_n;
    if (mu_n < mumin) mumin = mu_n;
    _contact_potential[it->second] = mu_n;

    id++;
  }

  double Ec = _ext_module->get_band_edge("Ec");
  double Ev = _ext_module->get_band_edge("Ev");
  params.ec = Ec;
  params.ev = Ev;

  params.deltaec = opt.DEc;
  params.deltaev = opt.DEv;

  params.emin = sol_opt.get_option("Emin", -mumax-opt.n_kT*kbT);
  params.emax = sol_opt.get_option("Emax", -mumin+opt.n_kT*kbT);
  params.estep = opt.Estep;

  // contour integration parameters
  for (unsigned int i = 0; i < 2; ++i)
  {
    params.np_n[i] = opt.Np_n[i];
    params.np_p[i] = opt.Np_n[i];
  }
  params.n_poles = opt.n_poles;

  // real axis integration parameters
  int Np_real = 0;
  Np_real = sol_opt.get_option("Np_real", Np_real);
  if (Np_real < 0)
  {
    Np_real = (params.emax - params.emin) / opt.deltaE;
    ostringstream os;
    os << "Setting real line integration points to " << Np_real;
    Messages::info(os.str());
  }
  params.np_real[0] = Np_real;

  // probably it does not make sense to have both?
  int contour_points = params.n_poles + opt.Np_n[0] + opt.Np_n[1];
  if ((params.np_real[0] > 0) && (contour_points > 0))
    Messages::warning("Are you sure you want both contour and real axis integration?");

  // 0 -> min(mu)
  // 1 -> max(mu)
  // x -> use all contacts (for real axis integration)
  params.min_or_max = 1;
  // TODO this is set internally in libnegf
  //if (Np_real > 0)
  //  params.min_or_max = 2;

  params.n_kt = opt.n_kT;




  // now hand parameters to library
  _libnegf->set_parameters(params);

  _libnegf->set_iteration(1);

  // set not to compute Device-Contact blocks
  _libnegf->device_contact_dm(0);

  // setup the structure
  /* calls negf_init_structure(handler, ncont,
   *                surfstart, surfend, contend, npl, plend, cblk)
   *
   * ncont = number of contacts
   * surfstart = indices of start of surface blocks
   * surfend  = indices of end of surface blocks (for first contact last device index)
   * contend = indices of last contact index
   * npl = number of principal layers
   * plend = last index of each PL
   */
  int n_vars = _ext_module->get_number_of_bands();
  _contend.resize(_qc_n_dofs.size());
  for (unsigned int i = 0; i < _contend.size(); ++i)
    _contend[i] = _qc_n_dofs[i]*n_vars;

  _surfstart.resize(_quantum_contacts.size());
  _surfstart[0] = _device_n_dofs*n_vars + 1;
  for (unsigned int i = 1; i < _surfstart.size(); ++i)
    _surfstart[i] = _contend[i-1] + 1;

  _surfend = _surfstart;
  for (auto&& s : _surfend)
    s -= 1;

  /*
  cerr << "surfstart = ";
  for (auto&& s : surfstart)
    cerr << s << " ";
  cerr << endl;

  cerr << "surfend = ";
  for (auto&& s : surfend)
    cerr << s << " ";
  cerr << endl;

  cerr << "contend = ";
  for (auto&& s : contend)
    cerr << s << " ";
  cerr << endl;
  */



  if (plot_solution("LDOS"))
  {
    _libnegf->init_ldos(_device_n_dofs*n_vars);
    // this will have libnegf give the whole LDOS back,
    // without any projection
    _libnegf->set_ldos_indices(0, vector<int>(1, -1));
  }

  //if (opt.verbosity > 60) _libnegf->partition_info();

}

void
Negf::finalize(void)
{
  // restore the original atomistic structure
  if ( _hamil_type == "etb" )
  {
    _ext_module->get_atomistic_structure()->reorder(_inv_perm);
  }
}

void
Negf::compute_current(void)
{
  _libnegf->set_verbose(opt.verbosity);

  //std::cout << _libnegf->current() << std::endl;

  current.clear();
  current.resize(2,0.0);

  //TODO: elCurrent. hlCurrent
  current[0] = _libnegf->current("eV","A");
  current[1] = current[0];
  cout << "curr from negf : " << current[0] << endl;
}

void
Negf::do_solve(void)
{

  setup_negf();

  if (plot_solution(elDensity) || plot_solution(hlDensity) ||
      plot_solution("LDOS"))
  {
    (_qdens_sys->solution)->zero();
  }

  if ( plot_solution(elDensity) )
  {
    Messages::info("Computing Electronic Density");

    if (get_options().has_submodel("k_integration_density"))
    {
      _which_integration = INTDENSITYEL;

      _k_int_density->solve();

      transfer_density(_k_int_density->get_solution(), "el");

    }
    else
    {
      _k_vec.zero();
      
      setup_hamil();

      unsigned int n_vars = _sys_H->n_vars();
      std::vector<double> density(_device_n_dofs * n_vars, 0.0);

      _libnegf->density(density, "el");

      transfer_density(density, "el");

      // write qdens on eldensity
      //NumericVector<Number>& qdens = *_qdens_sys->solution;
      //unsigned int p_id = _qdens_sys->variable_number("edens");
      //unsigned int p_dofs = _qc_n_dofs[_quantum_contacts.size()-1];
      //unsigned int start = p_id * p_dofs;
      //_eldensity.reserve(_device_n_dofs);
      //for (unsigned int i=0; i < _device_n_dofs; i++)
      //  _eldensity.push_back( qdens(start+i) );


      finalize();
    }
      
    //for (unsigned int i=0; i < _device_n_dofs; i++)
    //  std::cout<<_eldensity[i]<<std::endl;

    Messages::info("Density done");

  }

  if ( plot_solution("hlDensity") )
  {
    Messages::info("Computing Hole Density");

    if (get_options().has_submodel("k_integration_density"))
    {
      _which_integration = INTDENSITYHL;

      _k_int_density->solve();

      transfer_density(_k_int_density->get_solution(), "hl");

    }
    else
    {
      _k_vec.zero();
      
      setup_hamil();

      unsigned int n_vars = _sys_H->n_vars();
      std::vector<double> density(_device_n_dofs * n_vars, 0.0);

      _libnegf->density(density, "hl");

      transfer_density(density, "hl");

      // write qdens on hldensity
      //NumericVector<Number>& qdens = *_qdens_sys->solution;
      //unsigned int p_id = _qdens_sys->variable_number("hdens");
      //unsigned int p_dofs = _qc_n_dofs[_quantum_contacts.size()-1];
      //unsigned int start = p_id * p_dofs;
      //_hldensity.reserve(_device_n_dofs);
      //for (unsigned int i=0; i < _device_n_dofs; i++)
      //  _hldensity.push_back( qdens(start+i) );

      finalize();
    }
      
    //for (unsigned int i=0; i < _device_n_dofs; i++)
    //  std::cout<<_hldensity[i]<<std::endl;

    Messages::info("Density done");

  }


  if ( plot_solution("Current") )
  {
    Messages::info("Computing Current");

    if (get_options().has_submodel("k_integration_current"))
    {
      _which_integration = INTCURRENT;

      _k_int_current->solve();

      current = _k_int_current->get_solution();

    }
    else
    {
       _k_vec.zero();

       setup_hamil();

       compute_current();

       finalize();
    }

    double u = get_mesh_units();
    unsigned int dim = get_mesh().mesh_dimension();
    double area_factor;

    switch (dim)
    {
      case 1:
        area_factor = 1e-4/(u*u);  // k is in 1/u  
        break;                    // u = 1e-9 m => u = 1e-9 * 10^2 cm
      case 2:                     // 1/m = 1e-2/u 1/cm
        area_factor = 1e-2/u;
        break;
      case 3:
        area_factor = 1.0;
    }

    //get degeneracy of first band
    //const MeshBase& mesh = get_mesh();
    //MeshBase::const_element_iterator el = mesh.active_elements_begin();
    //const Elem* elem = *el;
    //NegfModel* negfmod = get_bulk_model<NegfModel>(elem);
    double deg = _ext_module->get_degeneracy();

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
      //std::cout<<"(negf) I= "<<current[id]<<std::endl;
      //std::cout<<"(negf) sgn, deg area "<<sign<<" "<<deg<<" "<<area_factor<<std::endl;

      _contact_current[it->second] = sign * deg * current[id] * area_factor;
      id++;
    }

    Messages::info("Current done");
    plot_globaldata();

  }
}



void
Negf::plot_LDOS(const std::vector<double>& energies,
    const std::vector<std::vector<double>>& ldos,
    const string& name_suffix)
{
  if (ldos.empty())
    return;

  string mod(name_suffix);
  if (!mod.empty()) mod = "_" + mod;

  string file = get_output_directory() + "/" + get_output_filename_prefix()
      + "_LDOS" + mod + TiberCad::get_filename_suffix() + ".m";
  ofstream of(file);

  size_t esteps = energies.size();

  unsigned int bands = _ext_module->get_number_of_bands();
  int npoints = _device_n_dofs * bands;

  of << "energy = [";
  for (auto&& erg : energies) of << erg << " ";
  of << "];\n";


  DofMap& dof_map = _sys_H->get_dof_map();
  std::vector<unsigned int> dof_indices;

  //set<double> coordinates;
  vector<double> coordinates(_device_n_dofs);

  // TODO what to do if mesh is distributed?
  MeshBase::const_element_iterator       nd     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator nd_end = this->active_local_elements_end();
  for ( ; nd != nd_end; ++nd)
  {
    const Elem* elem = *nd;
    dof_map.dof_indices(elem, dof_indices, 0);

    for (unsigned int n = 0; n < elem->n_nodes(); ++n)
    {
      unsigned int id = _inv_perm[dof_indices[n]];
      coordinates[id / bands] = elem->point(n)(0);
    }
  }


  of << "x = [";
  for (size_t i = 0; i < _device_n_dofs; ++i)
    of << coordinates[i] << " ";
  of << "];\n";

  of << "LDOS = [";
  for (int i = 0; i < esteps; i++)
  {
    for (int j = 0; j < npoints; j += bands)
    {
      double node_dos = 0.0;
      for (int b = 0; b < bands; ++b)
        node_dos += ldos[i][j + b];
      of << node_dos << " ";
    }
    of << "\n";
  }
  of << "];\n";

  //of << "x=1:" << _device_n_dofs << ";\n";
  //of << "pcolor(x, energy, log(abs(LDOS))), shading flat\n";
  of << "pcolor(x, energy, abs(LDOS)), shading flat\n";
  of << "ylabel('Energy')\n";
  of << "xlabel('x')\n";
}




void
Negf::occupy_LDOS(const std::vector<double>& ldos)
{
  // the following is quite ineffective, but it was easier like that...

  SimulationInterface* potmodel = SimulationInterface::find_simulation(opt.pot_module);
  ID efermi_id = INVALID_ID, hfermi_id = INVALID_ID;

  if (potmodel != NULL)
  {
    efermi_id = potmodel->get_solution_id("eQFermi");
    hfermi_id = potmodel->get_solution_id("hQFermi");
  }

  if ((efermi_id == INVALID_ID) || (hfermi_id == INVALID_ID))
  {
    Messages::warning("cannot occupy LDOS with classical Fermi levels");
    return;
  }


  if (!potmodel->is_solved())
      throw SolveFailedException("Simulation " + opt.pot_module +
          " must be solved first");


  double Ec = get_band_edge("Ec");
  double Ev = get_band_edge("Ev");

  // energy grid
  int esteps = (opt.Emax - opt.Emin) / opt.Estep + 1;
  unsigned int bands = _ext_module->get_number_of_bands();
  int npoints = _device_n_dofs * bands;


  int limit_id = (0.5 * (Ec + Ev) - opt.Emin) / opt.Estep;

  double u = get_mesh_units();
  unsigned int dim = get_mesh().mesh_dimension();

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
      break;
  }

  DofMap& dof_map_qdens = _qdens_sys->get_dof_map();
  std::vector<unsigned int> dof_indices_eldens;
  std::vector<unsigned int> dof_indices_hldens;

  NumericVector<Number>& qdens = *_qdens_sys->solution;
  qdens.zero();

  int el_id = _qdens_sys->variable_number("edens");
  int hl_id = _qdens_sys->variable_number("hdens");

  // compute connectivity of each node
  // we need the connectivity of the nodes to not double count
  //std::vector<int> connectivity(_device_n_dofs * n_vars, 0);
  //std::cout<<"qdens.size: "<<qdens.size()<<std::endl;

  std::vector<int> connectivity(qdens.size(), 0);
  {
    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_eldens, el_id);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        connectivity[dof_indices_eldens[n]]++;
    }
  }


  unsigned int n_vars = _sys_H->n_vars();

  vector<double> occupied_states(ldos.size());


  // The _sys_H system contains the nodal dofmap (all variables)
  DofMap& dof_map = _sys_H->get_dof_map();
  std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator el = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  NegfModel* negfmod = get_bulk_model<NegfModel>(*el);
  double deg = _ext_module->get_degeneracy();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map_qdens.dof_indices(elem, dof_indices_eldens, el_id);
    dof_map_qdens.dof_indices(elem, dof_indices_hldens, hl_id);

    // get the Fermi levels on the centroid
    double e_qfermi = 0;
    double h_qfermi = 0;
    potmodel->get_solution(elem, efermi_id, e_qfermi, elem->centroid());
    potmodel->get_solution(elem, hfermi_id, h_qfermi, elem->centroid());

    for (int band = 0; band < n_vars; band++)
    {
      dof_map.dof_indices(elem, dof_indices, band);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        // the index in the ldos vector
        unsigned int id = _inv_perm[dof_indices[n]];

        double edens = 0;
        double hdens = 0;
        for (int i = 0; i < esteps; i++)
        {
          if (i > limit_id)
          {
            double f = fermi(opt.Emin + i * opt.Estep, e_qfermi, SimulationOptions::T);
            edens += ldos[i + id*esteps] * f * opt.Estep;
            occupied_states[i + id*esteps] += ldos[i + id*esteps] * f / connectivity[dof_indices_eldens[n]];
          }
          else
          {
            double f = fermi(-(opt.Emin + i * opt.Estep), -h_qfermi, SimulationOptions::T);
            hdens += ldos[i + id*esteps] * f * opt.Estep;
            occupied_states[i + id*esteps] += ldos[i + id*esteps] * f / connectivity[dof_indices_eldens[n]];
          }
        }


        //std::cout<<band*n_dofs+n<<" dens: "<<dens<<std::endl;

        double vale = equ * deg * edens / connectivity[dof_indices_eldens[n]];
        double valh = equ * deg * hdens / connectivity[dof_indices_eldens[n]];

        qdens.add(dof_indices_eldens[n], vale);
        qdens.add(dof_indices_hldens[n], valh);
      }

    }
  }

  qdens.close();

  //plot_LDOS(occupied_states, "OCC");

}




void
Negf::calculate_for_k_point(const Point& k_point,
                                  DofField& field,
                                  double& error)
{

   for(short i=0;i<3;i++) _k_vec(i) = k_point(i);

   setup_hamil();

   //unsigned int n_vars = _sys_H->n_vars();
   //field.resize(_device_n_dofs * n_vars);
   unsigned int n_vars = _ext_module->get_H_dim();
   field.resize(n_vars);

   switch (_which_integration)
   {
     case INTDENSITYEL:
     {

       _libnegf->density(field, "el");

       //std::string out_file = "density_new.dat";
       //std::fstream ff(out_file.c_str(),std::fstream::out);
       //for (unsigned int i = 0; i < field.size(); i++)
       //      ff<< field[i] <<std::endl;
       //ff.close();

       error = 0.0;

       for (unsigned int i=0; i < field.size(); i++)
       {
         error += field[i];
       }

       error /= _device_n_dofs;

       //cout<<"(negf) density error: "<<error<<endl;

       break;
     }

     case INTDENSITYHL:
     {

       _libnegf->density(field, "hl");

       error = 0.0;

       for (unsigned int i=0; i < field.size(); i++)
       {
         error += field[i];
       }

       error /= _device_n_dofs;

       break;
     }


     case INTCURRENT:
     {
       field.clear();
       compute_current();

       // get_energies
       vector<double> erg;
       _libnegf->get_energies(erg);

       vector<vector<double>> transmission;
       // get transmission
       _libnegf->get_transmission(transmission);

       ostringstream os;
       os << "transmission";
       //if (k_point.norm() > 1e-6)
       os << "_k=(" << k_point(0) << ","
           << k_point(1) << "," << k_point(2) << ")";
       os << ".dat";

       ostringstream header;
       header << "# Transmission at k = "
              << "(" << k_point(0) << "," << k_point(1) << ","
              << k_point(2) << ")\n"
              << "#\n"
              << "# energy T\n";

       print_energy_resolved(os.str(), erg, transmission, header.str());

       //double curr = _libnegf->current();
       //_contact_potential[]
       //field.push_back(curr);
       //field.push_back(-curr);
       field = current;

       error = current[0];
       break;
     }
   }


   if (plot_solution("LDOS"))
   {
     Messages::info("Plot LDOS");

     // get_energies
     vector<double> erg;
     _libnegf->get_energies(erg);

     vector<vector<double>> ldos;
     _libnegf->get_ldos(ldos);

     ostringstream os;
     os << "k=(" << k_point(0) << ","
        << k_point(1) << "," << k_point(2) << ")";

     plot_LDOS(erg, ldos, os.str());
   }
      
   finalize();

}


void
Negf::print_energy_resolved(const string& file, const vector<double>& energy,
    const vector<vector<double>>& data, const string& header) const
{
  if (get_communicator().rank() != 0)
    return;

  ofstream out(get_output_directory() + "/" +
      get_output_filename_prefix() + "_" + file);

  if (!out.good())
      throw std::runtime_error("Could not open " + file + " for writing.");

  out << header;

  size_t N = energy.size();

  // sanity check on data
  for (size_t j = 0; j < data.size(); j++)
  {
    if (data[j].size() != N)
      throw std::runtime_error("Data is inconsistent with energy vector.");
  }

  for (size_t i = 0; i < N; ++i)
  {
    out << energy[i];
    for (size_t j = 0; j < data.size(); j++)
      out << " " << data[j][i];
    out << endl;
  }
}


void
Negf::transfer_density(const std::vector<double>& density, const std::string& particle)
{

  // compute total number of n_vars n_dofs
  unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_dofs = _qc_n_dofs[_quantum_contacts.size()-1];

  int particle_id;

  if (particle == "el") particle_id = _qdens_sys->variable_number("edens");
  if (particle == "hl") particle_id = _qdens_sys->variable_number("hdens");


  double equ;
  double u = get_mesh_units();
  unsigned int dim = get_mesh().mesh_dimension();
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
      break;
  }

  DofMap& dof_map_qdens = _qdens_sys->get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  // setup output vector qdens

  NumericVector<Number>& qdens = *_qdens_sys->solution;
  //qdens.zero();

  // compute connectivity of each node
  // we need the connectivity of the nodes to not double count
  //std::vector<int> connectivity(_device_n_dofs * n_vars, 0);
  //std::cout<<"qdens.size: "<<qdens.size()<<std::endl;

  std::vector<int> connectivity(qdens.size(), 0);
  {
    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
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

  MeshBase::const_element_iterator el = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  NegfModel* negfmod = get_bulk_model<NegfModel>(*el);

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map_qdens.dof_indices(elem, dof_indices_qdens, particle_id);

    for (int band = 0; band < n_vars; band++)
    {
      dof_map.dof_indices(elem, dof_indices, band);


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        //double dens = abs(density[_inv_perm[dof_indices[n]]]);
        double dens = density[_inv_perm[dof_indices[n]]];

        //std::cout<<band*n_dofs+n<<" dens: "<<dens<<std::endl;

        double val = equ * dens / connectivity[dof_indices_qdens[n]];

        qdens.add( dof_indices_qdens[n], val);
      }

      //for (unsigned int n = 0; n < elem->n_nodes(); n++)
      //std::cout<<connectivity[ dof_indices_qdens[n]]<<std::endl;

    }
  }

  qdens.close();
}


/*
void
Negf::calculate_density(const std::string& particle)
{
  double u = get_mesh_units();

  if (_ext_module != NULL)
  {
    const Scaling& sc = _ext_module->get_scaling();
    //u = sc.get_length_scaling() / sc.get_calc_mesh_units();
  }

  unsigned int dim = get_mesh().mesh_dimension();


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
 
  _libnegf->density(density, particle);

  // ---- debug ---------------------------------------------
  std::string out_file = "density_new.dat";
  std::fstream ff(out_file.c_str(),std::fstream::out);
  ff << "nvars: "<<n_vars <<endl;
  ff << "device ndofs: "<<_device_n_dofs <<endl;
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

  libMesh::DofMap& dof_map_qdens = _qdens_sys->get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  // setup output vector qdens
  libMesh::NumericVector<Number>& qdens = *_qdens_sys->solution;

  //qdens.zero();

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
  libMesh::DofMap& dof_map = _sys_H->get_dof_map();
  std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator el = get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();

  NegfModel* negfmod = get_bulk_model<NegfModel>(*el);
  double deg = _ext_module->get_degeneracy();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map_qdens.dof_indices(elem, dof_indices_qdens, particle_id);

    for (int band = 0; band < n_vars; band++)
    {
      dof_map.dof_indices(elem, dof_indices, band);


      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        double dens = abs(density[_inv_perm[dof_indices[n]]]);

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
*/


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

  opt.deltaE =  sol_opt.get_option("deltaE",opt.delta/2.0);

  int Np = 2 * opt.n_kT/opt.deltaE;

  opt.Np_real = sol_opt.get_option("Np_real",Np);

  opt.DEc = sol_opt.get_option("deltaEc", 0.25);

  opt.DEv = sol_opt.get_option("deltaEv", 0.25);

  opt.n_poles = sol_opt.get_option("Npoles", 0);

  opt.verbosity = get_option("verbosity",0);
  opt.verbosity = sol_opt.get_option("verbosity",0);
  //sol_opt.check_unused();
  opt.writeLDOS = plot_solution("LDOS");
  opt.writeLDOS = sol_opt.get_option("writeLDOS", opt.writeLDOS);

}

void
Negf::do_setup_solution_variables(void)
{


  declare_solution(ReorderPotential, REAL, NODES, "");
  declare_solution(elDensity, REAL, NODES, "1/cm^3");
  declare_solution(hlDensity, REAL, NODES, "1/cm^3");
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





bool
Negf::is_generalized(void)
{
  return _ext_module->is_generalized();
}

/*
 *'outputd10nml50/Hr.m'
 *'outputd10nml50/Hi.m'
 *'outputd10nml50/Sr.m'
 *'outputd10nml50/Si.m'
 *2                 ! n_cont
 *0                 ! n_blocks (0 libnegf do it)
 *26118 27270       ! cont_end 
 *24966 26118       ! surf_end
 *-2 0              ! Ec Ev
 *0 0               ! DEc DEv
 *-0.25 0 0.01      ! Emin Emax Estep
 *0.025 0.025       ! kT
 *1                 ! weight of k-point
 *20 20             ! Np_n(1:2)
 *0 0               ! Np_p(1:2)
 *10                ! Real axis points
 *3                 ! n_kT
 *3                 ! n_poles
 *2                 ! degeneracy factor (spin)
 *1e-05             ! delta
 *0                 ! n_ldos
 *                  ! --
 *0.115355 0.184647 ! mu_n  
 *0.115355 0.184647 ! mu_p 
 */

void
Negf::print_Lib()
{
  ModelOptions& sol_opt = get_solver_options();
  double DeltaEc = -opt.DEc;
  double DeltaEv =  opt.DEv;

  // Takes kb in eV/K so assuming eV as energy units.
  double kbT = SimulationOptions::temperature * Constants::kb;
  double wght = 1.0;


  // for now, use the same parameters
  std::vector <double> Np_p(2);
  for (unsigned int i = 0; i < 2; i++)
  {
    Np_p[i] = opt.Np_n[i];
  }

  // spin degeneracy must be passed 
  unsigned int n_vars = _ext_module->get_number_of_bands();
  double Ec = _ext_module->get_band_edge("Ec");
  double Ev = _ext_module->get_band_edge("Ev");

  unsigned int spin = _ext_module->get_degeneracy();
  unsigned int nLDOS = 0;
  unsigned int n_bands = _ext_module->get_number_of_bands();
  if (opt.writeLDOS) nLDOS = _device_n_dofs * n_bands;

  // NOTE: adapted to new libnegf working with indeces 
  std::vector <unsigned int> LDOSindeces(nLDOS);

  unsigned int ctr = 1;
  for (unsigned int i = 0; i < nLDOS; ctr++)
  {
    LDOSindeces[i] = ctr;
  }

  // phi: potential at boundaries (quantum contacts)
  // mu : electrochemical potential at boundaries (qc)
  std::vector <double> phi(_quantum_contacts.size(), 0.0);
  std::vector <double> mu_n(_quantum_contacts.size(), 0.0);
  std::vector <double> mu_p(_quantum_contacts.size(), 0.0);

  mumin=10000;
  mumax=-10000;
  ID id = 0;
  std::map<ID, QuantumContact*>::iterator it = _quantum_contacts.begin();
  const std::map<ID, QuantumContact*>::iterator end = _quantum_contacts.end();
  for (; it != end; ++it)
  {
    get_boundary_potentials(it->second, phi[id], mu_n[id], mu_p[id]);
     
    if(mu_n[id]>mumax){mumax = mu_n[id];}
    if(mu_n[id]<mumin){mumin = mu_n[id];}
    _contact_potential[it->second] = mu_n[id];

    id++;
  }

  opt.Emin = sol_opt.get_option("Emin",-mumax-opt.n_kT*kbT);
  opt.Emax = sol_opt.get_option("Emax",-mumin+opt.n_kT*kbT);

  std::string outpath = get_scratch_directory();
  std::string out_file = outpath + "/negf.in";
  std::fstream ff(out_file.c_str(),std::fstream::out);

  // we pass matrices by memory
  // TODO libnegf reads two strings, first is discarded: what is it for?
  ff << "Hr 'memory'" << endl;
  ff << "Hi 'memory'" << endl;
  ff << "Sr 'memory'" << endl;
  ff << "Si 'memory'" << endl;
  /*
  ff<<"'"+outpath+"/Hr.m'"<<std::endl;
  ff<<"'"+outpath+"/Hi.m'"<<std::endl;
  
  if (is_generalized())
  {  
    ff<<"'"+outpath+"/Sr.m'"<<std::endl;
    ff<<"'"+outpath+"/Si.m'"<<std::endl; 
  }
  else
  {
    ff<<"'identity'"<<std::endl;
    ff<<"'identity'"<<std::endl;
  }
  */
  
  ff<<"N.conts  "<<_quantum_contacts.size()<<std::endl;
  
  
  // Nblocks (if 0 are partitioned by the library)
  unsigned int nPL=sol_opt.get_option("number_of_PL",0);
  ff<<"N.PLs  "<<nPL<<std::endl;
  if (nPL>0)
  { 
    unsigned int sizePL=sol_opt.get_option("size_of_PL",1);
    for (unsigned int i=0; i<nPL; i++)
       ff<<(i+1)*sizePL*n_vars<<" ";

    ff<<std::endl;
  }

  //ff<<floor(_device_n_dofs*n_vars/2)<<" "<<_device_n_dofs*n_vars<<std::endl;

  // contact end dofs:
  ff<<"Cont.dofs ";
  for (unsigned int i = 0; i <_quantum_contacts.size(); i++)
    ff<<_qc_n_dofs[i]*n_vars<<" ";
  ff<<std::endl;

  ff<<"Dev.dofs  "<<_device_n_dofs*n_vars<<" ";
  for (unsigned int i = 0; i <(_quantum_contacts.size()-1); i++)
    ff<<_qc_n_dofs[i]*n_vars<<" ";
  ff<<std::endl;

  ff<<"Emin.Emax  "<<opt.Emin<<" "<<opt.Emax<<" "<<opt.Estep<<std::endl;
  ff<<"kbT  "<<kbT<<"  "<<kbT<<std::endl;
  ff<<"wght  "<<wght<<std::endl;

  ff<<"Np_n  ";
  for (unsigned int i = 0; i < 2; i++)
    ff<<opt.Np_n[i]<<" ";
  ff<<std::endl;

  ff<<"Np_p  ";
  for (unsigned int i = 0; i < 2; i++)
    ff<<Np_p[i]<<" ";
  ff<<std::endl;

  int Np;
  if (sol_opt.find_option("Np_real"))
  {
     Np = opt.Np_real;
  }
  else
  {  
     Np = (abs(mu_n[0]-mu_n[1])+ 2 * opt.n_kT)/opt.deltaE ;
     if (Np>1000) 
     {	   
        std::ostringstream os;
        os << "Np_real has been set to "<<Np;
	      Messages::warning(os.str());
     }
  }

  ff<<"Np_real  "<< Np <<std::endl;
  ff<<"n_kT  "<<opt.n_kT<<std::endl;
  ff<<"n_poles  "<<opt.n_poles<<std::endl;
  ff<<"spin  "<<spin<<std::endl;
  ff<<"delta  "<<opt.delta<<std::endl;
  ff<<"nLDOS  "<<nLDOS<<std::endl;
  
  if (nLDOS>0)
  {
    ff<<"LDOS  ";
    for (unsigned int i = 0; i < nLDOS; i++)
      ff<<LDOSindeces[i]<<" ";

    ff<<std::endl;
  }

  ff<<"Ec.Ev  ";
  
  for (unsigned int i = 0; i < nLDOS; i++)
    ff<<Ec<<" "<<Ev<<std::endl;
  
  ff<<std::endl;
  // takes mu_n (ELECTROCHEMICAL potential for electrons) 
  //       mu_p (ELECTROCHEMICAL potential for holes)
  // NOTE: Ef has been removed now. 
  // libNEGF works with eletrochemical potential, mu
  // (e=|e|) For Ef=0, mu_n=-eV, mu_p=eV:
  // OLD: Ef-(-mu_n) = mu_n ; Ef-mu_p = -mu_p
  // NEW: pass mu_n and -mu_p
  ff<<"mu_n  ";
  for (unsigned int i = 0; i < _quantum_contacts.size(); i++)
    ff<<mu_n[i]<<" ";
  ff<<std::endl;

  ff<<"mu_p  ";
  for (unsigned int i = 0; i < _quantum_contacts.size(); i++)
    ff<<-mu_p[i]<<" ";
  ff<<std::endl;
  ff<<"DEc.DEv  "<<DeltaEc<<" "<<DeltaEv<<std::endl;    
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
    const libMesh::NumericVector<Number>& solution = _sys->get_solution_vector();
    const libMesh::DofMap& dof_map = _sys->get_dof_map();
    ID u_var = _sys->variable_number("u0");
    FEType fe_type = _sys->variable_type(u_var);
    UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));
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

  bool do_edens = values.count(elDensity);
  bool do_hdens = values.count(hlDensity);
  if (do_edens || do_hdens)
  {
    //std::cout<<"(Negf) Getting electron density"<<std::endl;
    _sys_H = &get_equation_system<TiberLinearSystem>(0);
    // internal variables of _qdens_sys (el: 0, hl: 1)
    ID el_id = _qdens_sys->variable_number("edens");
    ID hl_id = _qdens_sys->variable_number("hdens");

    libMesh::NumericVector<Number>& qdens = *_qdens_sys->solution;

    FEType fe_type = _qdens_sys->variable_type(el_id);
    UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &p);

    libMesh::DofMap& dof_map = _qdens_sys->get_dof_map();
    std::vector<unsigned int> dof_indices_e;
    std::vector<unsigned int> dof_indices_h;
    dof_map.dof_indices(elem, dof_indices_e, el_id);
    dof_map.dof_indices(elem, dof_indices_h, hl_id);
    unsigned int n_dofs = phi.size();

    // the phi^2 factor comes from the fact that the more correct interpolation is the square
    // of the basis function, because the probability densities are the square of the states
    // NOTE: maybe one should check if this gives really a better result
    // NOTE: 2011-12-01 the above turned out to be wrong
    //  phi[i] rho[i][j] phi[j] ~~  phi[i] rho[i][i] phi[i]
    // rho[i][i] is the correct density at node i: phi is used for linear interpolation.

    if (do_edens)
    {
      for (unsigned int n = 0; n < np; n++)
      {
        double value_n = 0;
 
        for (unsigned int i = 0; i < n_dofs; i++)
          value_n += phi[i][n] * qdens(dof_indices_e[i]);

        values[elDensity][n] = value_n;

      }
    }

    if (do_hdens)
    {
      for (unsigned int n = 0; n < np; n++)
      {
        double value_p = 0;

        for (unsigned int i = 0; i < n_dofs; i++)
          value_p += phi[i][n] * qdens(dof_indices_h[i]);

        values[hlDensity][n] = value_p;
      }
    }
  }

  if (values.count(eCurrentDensity))
  {
    //std::cout<<"(Negf) Getting e current density"<<std::endl;
    //given an element/point decide which contact it belongs to
    //values[eCurrentDensity] = std::vector<double>(np*3, 0.0);

    // loop over the sides and check on which boundary it is
    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      const ElementSide elside(elem, side);
      Boundary* b = _env->get_boundary(elside);
      QuantumContact* qc = _qc_boundaries[b];

      if (b != NULL && qc != NULL)
      {
        Point point_current;      
        if (dim>1)
        { 
          UniquePtr<Elem> elside = elem->build_side(side); 		
          point_current = _contact_current[qc] * qc->get_normal() * elside->volume() / qc->get_area();
        }
        else
        {
          point_current = _contact_current[qc] * qc->get_normal(); 
        }
	      
        // assign the same current to all points
        for (unsigned int n = 0; n < np; n++)
        {
	         //std::cout<<"point: "<<p[n](0)<< " "<<p[n](1)<<" "<<p[n](2)<<std::endl;
           values[eCurrentDensity][n*3]   =  point_current(0);
           values[eCurrentDensity][n*3+1] =  point_current(1);
           values[eCurrentDensity][n*3+2] =  point_current(2);
	         //std::cout<<"current: "<<point_current(0)<< " "<<point_current(1)<<" "<<point_current(2)<<std::endl;
        }
        break;
      }

    }
    //std::cout << "done" <<std::endl;
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


// -----------------------------------------------------------------------
// Assemble a Laplace problem in active regions and use solution to
// reorder dofs
void
Negf::reorder(void)
{

  //
  // The following does only make sense if the the mesh in this module is
  // the same as the one of the Hamiltonian provider
  //
  if (&this->get_mesh() != &_ext_module->get_mesh())
    throw InitFailedException("mesh in NEGF and Hamiltonian provider must be the same.");

  unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_var_reord = 1;
  //std::cerr<<"n vars "<<n_vars<<std::endl;

  std::vector<ID> u_vars(n_var_reord,0);
  set<ID> qc_ids;
  for (auto&& id : _quantum_contacts)
    qc_ids.insert(id.first);

  qc_ids.insert(get_region_ids().begin(), get_region_ids().end());

  // now qc_ids contains region ids of device and quantum contacts

  for (unsigned int n=0; n< n_var_reord; n++)
  {
    std::ostringstream var_str;
    var_str << "u" << n;
    std::string name = var_str.str();
    _sys->add_variable(name, FIRST, LAGRANGE, &qc_ids);
    u_vars[n] = _sys->variable_number(name);
  }

  _sys->attach_assemble_object(_reorder_assembly);

  _sys->init();

  _sys->solve();

  //std::cerr<<"Laplace solved"<<std::endl;

  const libMesh::NumericVector<Number>& solution = _sys->get_solution_vector();

  const libMesh::DofMap& dof_map = _sys->get_dof_map();

  std::vector<unsigned int> dof_indices_u;

  unsigned int sol_size = solution.size();
  //cerr << "solution size = " << sol_size << endl;

  // setup initial permutation vector as identitiy
  // the vector runs only on the device region where the
  // dof reordering is performed
  _perm.resize(_device_n_dofs, 0);
  for (unsigned int i = 0; i < _device_n_dofs; i++)
    _perm[i] = i;

  // Number of dofs of the first band
  unsigned int n_dofs = _qc_n_dofs[_quantum_contacts.size()-1];

  //std::cerr<<"Sorting"<<std::endl;
  //std::cerr<<"n_dofs: "<<n_dofs<<std::endl;
  //std::cerr<<"dev: "<<_device_n_dofs<<std::endl;


  //
  // NOTE: it seems there is a difference between > and < here,
  // maybe something related to the reference contact?
  class Compare
  {
    public:
      Compare(const libMesh::NumericVector<Number>& v) : _v(v) {};
      bool operator()(size_t i, size_t j)
      {
        return(_v(j) > _v(i));
      }

    private:
      const libMesh::NumericVector<Number>& _v;
  };

  std::sort(_perm.begin(), _perm.end(), Compare(_sys->get_solution_vector()));

  //
  // Now _perm contains the permutation of the nodes, basically, and we can now calculate
  // the PL blocks. For that we need to get the number of DOFs on each node of the external
  // module. At the same time, we prepare the DOF permutation tables.

  unsigned int nPL = get_solver_options().get_option("number_of_PL", 0);
  unsigned int size_of_PL = get_solver_options().get_option("PL_size", 16);

  ostringstream os;
  os << "Number of principal layers (PL) ";

  if (nPL < _quantum_contacts.size())
  {
    nPL = _quantum_contacts.size();
    os << "(automatic choice)";
  }
  os << ": " << nPL;
  Messages::info(os.str());

  _end_blocks.resize(0);
  _end_blocks.resize(nPL, 0);

  vector<double> ranges(nPL);
  for (unsigned int i = 0; i < nPL; ++i)
    ranges[nPL-1-i] = (_quantum_contacts.size() - 1) *
                           (1.0 - static_cast<double>(i)/nPL);


  unsigned int sysid = _ext_module->get_equation_system_id();
  unsigned int mysid = _sys->number();
  unsigned int n_dofs_total = 0;
  //cerr << "sysid : " << sysid << endl;
  //cerr << "_device_n_dofs : " << _device_n_dofs << endl;

  auto it = get_mesh().active_nodes_begin();
  const auto end = get_mesh().active_nodes_end();
  for ( ; it != end; ++it)
  {
    const Node* node = *it;
    if (node->n_dofs(mysid, 0))
    {
      unsigned int dof = node->dof_number(mysid, 0, 0);

      if (dof < _device_n_dofs)
      {
      //cerr << *node << endl;
        // get the number of dofs in the Hamiltonian
        unsigned int n_hamil_dofs = node->n_dofs(sysid);

        //unsigned int real_dof = _perm[dof];

        // the reorder "potential"
        double v = _sys->get_solution_vector()(dof);

        //find the position in the ranges and add number of dofs
        unsigned int i = 0;
        for ( ; (i < nPL) && (v > ranges[i]); ++i);
        //cerr << v << " -> " << i << " # dof " << n_hamil_dofs << endl;
        _end_blocks[i] += n_hamil_dofs;
        n_dofs_total += n_hamil_dofs;
      }
    }
  }

  for (unsigned int i = 1; i < nPL; ++i)
    _end_blocks[i] += _end_blocks[i-1];


  os.str("");
  os << "# DOFs       : " << n_dofs_total << "\n";
  os << "mean PL size : " << n_dofs_total/nPL << "\n";
  Messages::info(os.str());

  /*
  if (verbose() > 3)
  {
    cerr << "PLs (block end indices):\n";
    for (auto&& b : _end_blocks)
      cerr << b << " ";
      cerr << endl;
  }
  */


  // at last, prepare permutation

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

  // NOTE 2021-09-30 since libmesh 1.0 DOFs are organized in node-major order
  for (unsigned int i = 0; i < _device_n_dofs ; i++)
    for (unsigned int k = 0; k < n_vars; k++)
      _inv_perm[_perm[i]*n_vars+k]= i*n_vars+k;

  // Reset dofs in QC as an identity
  for (unsigned int i = 0; i < n_dofs - _device_n_dofs ; i++)
    for (unsigned int k = 0; k < n_vars; k++)
    {
      size_t index = (_device_n_dofs + i)*n_vars + k;
      _inv_perm[index] = index;
    }

  // we recalculate the full permutation vector
  _perm.resize(_inv_perm.size());
  for (unsigned int i = 0; i < _perm.size(); ++i)
    _perm[_inv_perm[i]] = i;

  //for (auto&& p : _perm)
  //  cerr << p+1 << " ";
  //cerr << endl;

  //for (auto&& p : _inv_perm)
  //  cerr << p+1 << " ";
  //cerr << endl;
}



//Approch to Laplace problem
void
Negf::reorder_assemble(void)
{

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

  libMesh::DofMap& dof_map = _sys->get_dof_map();

  unsigned int n_vars = _sys->n_vars();

  libMesh::FEType fe_type = dof_map.variable_type(0);

  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));

  libMesh::QGauss qrule(dim, THIRD);

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  libMesh::DenseMatrix<Number> Ke;
  libMesh::DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  const double penalty = 1.e10;

  for (unsigned int k=0; k < n_vars; k++)
  {
    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();


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
          ID cc = _bd_num[bd];

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

  // Iterates over QuantumContact regions by their subdomain ids
  set<ID> qc_ids;
  for (auto&& id : _quantum_contacts)
  {
    qc_ids.insert(id.first);
  }


  for (unsigned int k=0; k < n_vars; k++)
  {
    auto       el     = mesh.active_subdomains_elements_begin(qc_ids);
    const auto end_el = mesh.active_subdomains_elements_end(qc_ids);

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

        unsigned int cc = _bd_num[_env->get_boundary(qc_name)];  // cc = 0 or 1
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



void
Negf::get_boundary_potentials(QuantumContact* qc, double& av_V, double& av_mu_n, double& av_mu_p)
{
  PotentialInterface pot_model;
  PotentialInterface mue_model;
  PotentialInterface muh_model;
  av_V = 0.0;
  av_mu_n = 0.0;
  av_mu_p = 0.0;

  if (opt.pot_module == "none") return;
  else
  {
    pot_model.set_simulation(opt.pot_module+".ElPotential");
    mue_model.set_simulation(opt.pot_module+".eQFermi");
    muh_model.set_simulation(opt.pot_module+".hQFermi");
  }

  std::vector<double> V;

  MeshBase& mesh = get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  UniquePtr<FEBase> fe( FEBase::build(dim, FEType() ));

  libMesh::QGauss qrule(dim, FIRST); // Why order 1 ?

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

        av_V += pot_model.get_potential(pair.first, pair.second) * elem->volume()/volume ;
        av_mu_n += mue_model.get_potential(pair.first, pair.second) * elem->volume()/volume;
        av_mu_p += muh_model.get_potential(pair.first, pair.second) * elem->volume()/volume;
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




const Boundary*
Negf::get_boundary(const QuantumContact* qc)
{
  std::map<const QuantumContact*, const Boundary*>::iterator it = _bd_map.find(qc);
  return it->second;
}


