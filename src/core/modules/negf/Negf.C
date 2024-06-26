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
#include "Kspace.h"
#include <mpi.h>


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
  
  // reorder contacts according to minimum coordinate
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
  MeshBase::const_element_iterator el = active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = active_local_elements_end();
  NegfModel* negfmod;

  negfmod = get_bulk_model<NegfModel>(*el);

  if (negfmod == nullptr)
    throw InitFailedException("Cannot find bulk model in NEGF");

  _hamil_type = negfmod->get_model_name(0);

  std::cout<<"(negf) init: "<< _hamil_type  <<std::endl;

  std::string sim = negfmod->get_simulation(0);
  std::cout<<"(negf) sim: "<< sim <<std::endl;

  _ext_module = dynamic_cast<EigenvalueProblem*>(find_simulation(sim));

  if (_ext_module == nullptr)
  {
    throw InitFailedException("NEGF module needs an external"
        " provider of the Hamiltonian.");
  }
  else
  {
    if (!_ext_module->is_initialized())
      _ext_module->init();
  }

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

  // We add a second system just to contain the density
  if (plot_solution(elDensity) || plot_solution(hlDensity) ||
      plot_solution("LDOS"))
  {
    std::cout<<"(negf) create eq sys for elDensity"<<std::endl;
    ID id = create_equation_system("linear","");
    _qdens_sys = &get_equation_system<TiberLinearSystem>(id);
    _qdens_sys->add_variable("edens", libMeshEnums::FIRST, LAGRANGE, &get_region_ids());
    _qdens_sys->add_variable("hdens", libMeshEnums::FIRST, LAGRANGE, &get_region_ids());
    _qdens_sys->init();
  }
    
  
  std::cout<<"(negf) init done. " <<std::endl;
  
}


//! The initialization
void
Negf::init_efa_hamil(void)
{

  unsigned int n_bands;

  n_bands = _ext_module->get_number_of_bands();

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

   unsigned int nPL = get_solver_options().get_option("number_of_PL", 1);

   _end_blocks.resize(0);
   _end_blocks.resize(nPL, 0);

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

   //Order the atoms
   sortclass sortobj(atoms);
   std::sort(_perm.begin(),_perm.begin()+last_dev+1,sortobj);

   //for (unsigned int i=0; i< N_atoms; i++)
   //{
   //   std::cout<<"atom "<<i<<" reg "<<atoms[i].get_region_ID()
   //            <<" iperm "<<_inv_perm[i]<<" perm "<<_perm[i]<<std::endl;
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
   
   for (ID i=0; i<2; i++)
   {
      std::cout<<"Conts: "<<i+1<<"  "<<_qc_n_dofs[i]+1<<std::endl;
   }

   // Reorder the qc_n_dofs such that they are with increasing order
   for (unsigned int i = 0; i <(_quantum_contacts.size()-1); i++)
    if (_qc_n_dofs[i]>_qc_n_dofs[i+1])
      std::swap(_qc_n_dofs[i],_qc_n_dofs[i+1]);


   //Create the _end_blocks vector.
   unsigned int n_orbitals = _ext_module->get_number_of_bands(); 

   //_device_n_dofs is equal to N_atoms in the device
   unsigned int size_of_PL = _device_n_dofs / nPL;

   for (unsigned int i = 0; i <nPL; i++)
    _end_blocks[i] = size_of_PL * n_orbitals;
   
   unsigned int remaining_atoms = _device_n_dofs - nPL * size_of_PL;
   if (remaining_atoms != 0) _end_blocks[nPL-1] += remaining_atoms * n_orbitals;
   
   for (unsigned int i = 1; i < nPL; ++i)
    _end_blocks[i] += _end_blocks[i-1];

   // These two lines remove any duplicate (if they are for some reason created)
   auto bl_it = std::unique(_end_blocks.begin(), _end_blocks.end());
   _end_blocks.resize(distance(_end_blocks.begin(), bl_it));
   
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
           // Get lattice vectors for _init_basis routine
           if (_scattering) _lattice_vectors = get_lattice_vectors(r1, r2);
           break;

         case 1:
           r1[0] = cf*vectors[2](0); r1[1] = cf*vectors[2](1); r1[2] = cf*vectors[2](2);
           kopts.set_option("r1",r1);
           if (_scattering) _lattice_vectors = get_lattice_vectors(r1);
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
          if (_scattering) 
          {
            std::vector<double> r1(3);
            for (int i=0; i<3; i++) r1[i] = c(i);
            _lattice_vectors = get_lattice_vectors(r1);
          }
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
          if (_scattering) 
          {
            std::vector<double> r1(3);
            std::vector<double> r2(3);
            for (int i=0; i<3; i++)
            {
              r1[i] = b(i);
              r2[i] = c(i);
            } 
            _lattice_vectors = get_lattice_vectors(r1, r2);
          }
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

  // Override the k_dim derived from mesh dimension, if it was specified in input
  if (get_option("x-periodicity", false) || get_option("y-periodicity", false) || get_option("z-periodicity", false))
  {
    std::vector<bool> input_periodicity(3, false);
    input_periodicity[0] = get_option("x-periodicity", false);
    input_periodicity[1] = get_option("y-periodicity", false);
    input_periodicity[2] = get_option("z-periodicity", false);

    k_dim = std::accumulate(input_periodicity.begin(), input_periodicity.end(), 0);
  }

  // For ETB we can get the k-dimension directly from the atomistic structure (periodicity information)
  if ( _hamil_type == "etb" ) 
  {
    unsigned int k_dim_etb;
    std::vector<bool> periodicity = _ext_module->get_atomistic_structure()->get_periodicity_vector();
    k_dim_etb = std::accumulate(periodicity.begin(), periodicity.end(), 0);
    // Check if it was defined in the input file
    if (k_dim_etb != 0) k_dim = k_dim_etb;
  }

  k_dim = min(k_dim, 2u); // It cannot really be more than 2 with NEGF

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
                                              _k_comm);

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
                                               _k_comm);

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
    {
      _ext_module->get_H_csr(A, JA, IA);
    }

    else
      _ext_module->get_H_csr(A, JA, IA, _perm);

    if (_scattering)
    {
      _libnegf->create_HS_container(_n_Hk);
      _libnegf->set_H_csr(nrow, A, JA, IA, _iK+1); // +1 because Fortran is 1-indexed
    }
    else
    {
      _libnegf->create_HS_container(1);
      _libnegf->set_H_csr(nrow, A, JA, IA, 1);
    }
    

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


    if (_scattering)
    {
      _libnegf->set_S_csr(nrow, S, JS, IS, _iK+1);
    }
    else
    {
      _libnegf->set_S_csr(nrow, S, JS, IS, 1);
    }

    Messages::info("done.");
  }
  else
  {
    if (_scattering)
    {
      _libnegf->set_S_id(nrow, _iK+1);
    }
    else
    {
      _libnegf->set_S_id(nrow, 1);
    }
  }

  // we have to reinitialize some structures in libnegf
  if (get_option("print_matrices", false))
    _libnegf->print_mat();

  _cblk = _libnegf->contact_blocks(_quantum_contacts.size(),
                               _surfstart, _surfend, _contend,
                               _end_blocks.size(), _end_blocks);

  _libnegf->init_structure(_quantum_contacts.size(),
                           _surfstart, _surfend, _contend,
                           _end_blocks.size(), _end_blocks, _cblk);
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
  // first get the parameter structure with defaults
  NegfWrapper::Parameters params;
  _libnegf->get_parameters(params);

  // now set values
  params.verbose = opt.verbosity;
  //params.readOldDM_SGFs = ;
  //params.readOldT_SGFs = ;
  params.g_spin = _ext_module->get_degeneracy();
  {
    ostringstream os;
    os << "spin degeneracy : " << params.g_spin;
    Messages::info(os.str());
  }
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
    params.mu_p[id] =  mu_p;
    //params.mu[id] = ; // not needed, for DFTB
    // TODO a hack for now
    params.mu[id] = mu_n;
    
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
  params.estep_coarse = opt.Estep_coarse;

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
  params.np_real = Np_real;

  // probably it does not make sense to have both?
  int contour_points = params.n_poles + opt.Np_n[0] + opt.Np_n[1];
  if ((params.np_real > 0) && (contour_points > 0))
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
  /*
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
  for (auto&& s : _surfstart)
    cerr << s << " ";
  cerr << endl;

  cerr << "surfend = ";
  for (auto&& s : _surfend)
    cerr << s << " ";
  cerr << endl;

  cerr << "contend = ";
  for (auto&& s : _contend)
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

  if (_scattering)
  {
    // This works only after init_k_space_integration()
    init_basis();

    //Set orbsperatom in elastic scattering, only if some interactions were defined
    if (_interactions.size() > 0 && _interactions[0].model != DUMMY)
    {
      std::vector<int> orbitals_per_atom(_device_n_dofs, n_vars);

      for (int i = 0; i < _interactions.size(); i++)
      {
        Interaction* inter = &_interactions[i];
        if (inter->model == DEPHATOMBLOCK || inter->model == DEPHOVERLAP)
        {
          inter->orbsperatm.resize(_device_n_dofs, n_vars);
          inter->orbsperatm = orbitals_per_atom;
        }
      }
    }

    setup_interactions();
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

std::vector<double>
Negf::compute_layer_currents(void)
{
  _libnegf->set_verbose(opt.verbosity);

  int nPL = _end_blocks.size();
  std::vector<double> layer_current(nPL-1);

  _libnegf->layer_current(layer_current, "eV", "A");
  std::cout << "Layer currents from negf :" << std::endl;
  for (int i=0; i<layer_current.size(); i++)
  {
    std::cout << layer_current[i] << std::endl;
  }
  return layer_current;
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
    if (get_option("quasi_equilibrium",false))
      Messages::info("Quasi-equilibrium approximation");

    if (get_options().has_submodel("k_integration_density"))
    {
      if (_scattering)
      {
        unsigned int n_vars = _ext_module->get_number_of_bands();
        std::vector<double> density(_device_n_dofs * n_vars, 0.0);

        set_kpoints("elDensity");
        set_hamiltonians();
        // _libnegf->compute_density_inelastic(density, "el");
        // transfer_density();
      }
      else
      {
        _which_integration = INTDENSITYEL;
        _k_int_density->solve();
        transfer_density(_k_int_density->get_solution(), "el");
      }

    }
    else
    {
      if (_scattering)
      {
        Messages::warning("Scattering block was added without density k-integration.");
      }
      
      _k_vec.zero();

      setup_hamil();

      unsigned int n_vars = _ext_module->get_number_of_bands();
      std::vector<double> density(_device_n_dofs * n_vars, 0.0); //device_n_dofs = n nodi del device

      if (get_option("quasi_equilibrium",false))
      {
        std::vector<double> Ec(_device_n_dofs*n_vars);
        std::vector<double> Ev(_device_n_dofs*n_vars);
        std::vector<double> muN(_device_n_dofs*n_vars);
        std::vector<double> muP(_device_n_dofs*n_vars);

        get_mu_and_bands(Ec, Ev, muN, muP);

        _libnegf->quasi_equilibrium_density(density, "el", Ec, Ev, muN, muP);

      }
      else
      {
        _libnegf->density(density, "el");
      }
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
    if (get_option("quasi_equilibrium",false))
          Messages::info("Quasi-equilibrium approximation");

    if (get_options().has_submodel("k_integration_density"))
    {
      if (_scattering)
      {
        unsigned int n_vars = _ext_module->get_number_of_bands();
        std::vector<double> density(_device_n_dofs * n_vars, 0.0);

        set_kpoints("hlDensity");
        set_hamiltonians();
        // _libnegf->compute_density_inelastic(density, "hl");
      }
      else
      {
        _which_integration = INTDENSITYHL;
        _k_int_density->solve();
        transfer_density(_k_int_density->get_solution(), "hl");
      }
    }
    else
    {
      if (_scattering)
      {
        Messages::warning("Scattering block was added without density k-integration.");
      }

      _k_vec.zero();
      
      setup_hamil();

      // unsigned int n_vars = _sys_H->n_vars();
      unsigned int n_vars = _ext_module->get_number_of_bands();
      std::vector<double> density(_device_n_dofs * n_vars, 0.0);

      if (get_option("quasi_equilibrium",false))
      {
        std::vector<double> Ec(_device_n_dofs*n_vars);
        std::vector<double> Ev(_device_n_dofs*n_vars);
        std::vector<double> muN(_device_n_dofs*n_vars);
        std::vector<double> muP(_device_n_dofs*n_vars);

        get_mu_and_bands(Ec, Ev, muN, muP);

        _libnegf->quasi_equilibrium_density(density, "hl", Ec, Ev, muN, muP);

      }
      else
      {
        _libnegf->density(density, "hl");
      }
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
    if(_scattering)
    {
      Messages::info("Computing Layer Current");
      set_kpoints("Current");
      set_hamiltonians();

      std::vector<double> layer_curr = compute_layer_currents();

      if (plot_solution("LDOS"))
      {
        Messages::info("Plot LDOS");

        // get_energies
        vector<double> erg;
        _libnegf->get_energies(erg);

        vector<vector<double>> ldos;
        _libnegf->get_ldos(ldos);

        plot_LDOS(erg, ldos, "");
      }
      
      current.clear();
      current.resize(2,0.0);
      current[0] = layer_curr[0];
      current[1] = layer_curr[layer_curr.size()-1];
    }
    else
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
    }

    double area_factor = 1.0;

    // k is defined in units of 1/nm
    switch (get_mesh().mesh_dimension())
    {
      case 1:
        area_factor = 1e14;
        break;
      case 2:
        area_factor = 1e7;
        break;
      //case 3: // it remains 1 in this case
      default:
        break;
    }

    //get degeneracy of first band
    //const MeshBase& mesh = get_mesh();
    //MeshBase::const_element_iterator el = mesh.active_elements_begin();
    //const Elem* elem = *el;
    //NegfModel* negfmod = get_bulk_model<NegfModel>(elem);

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

      _contact_current[it->second] = sign * current[id] * area_factor;
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

  // MATLAB plotting commented, replaced with python
  string file = get_output_directory() + "/" + get_output_filename_prefix()
      // + "_LDOS" + mod + TiberCad::get_filename_suffix() + ".m";
      + "_LDOS" + mod + TiberCad::get_filename_suffix() + ".py";

  ofstream of(file);

  size_t esteps = energies.size();

  unsigned int bands = _ext_module->get_number_of_bands();
  int npoints = _device_n_dofs * bands;

  // of << "energy = [";
  // for (auto&& erg : energies) of << erg << " ";
  // of << "];\n";
  of << "import matplotlib.pyplot as plt\nimport numpy as np\n\n";
  of << "energy = np.array([";
  for (auto&& erg : energies) of << erg << ", ";
  of << "])\n";

  vector<double> coordinates(_device_n_dofs);

  // TODO what to do if mesh is distributed?
  if (_hamil_type == "efa")
  {
    DofMap& dof_map = _sys_H->get_dof_map();
    std::vector<unsigned int> dof_indices;

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
  }
  // Cannot use _sys_H for ETB, define coordinates using atomistic structure
  if (_hamil_type == "etb")
  {
    AtomisticStructure* as = _ext_module->get_atomistic_structure();
    unsigned int N_atoms = _device_n_dofs;
    std::vector<Atom>& atoms = as->get_structure_atoms();

    for (int i=0; i<N_atoms; i++)
    {
      // x coord from Angstrom to mesh units
      double equ = 1e-10 / get_mesh_units();
      coordinates[i] = atoms[i].get_position(0) * equ;
    }
  }

  //of << "x = [";
  of << "x = np.array([";

  for (size_t i = 0; i < _device_n_dofs; ++i)
    of << coordinates[i] << ", ";
  of << "])\n";
  //   of << coordinates[i] << " ";
  // of << "];\n";

  //of << "LDOS = [";
  of << "LDOS = np.array([\n";

  for (int i = 0; i < esteps; i++)
  {
    of << "["; // <- only for python plotting
    for (int j = 0; j < npoints; j += bands)
    {
      double node_dos = 0.0;
      for (int b = 0; b < bands; ++b)
        node_dos += ldos[i][j + b];
      of << node_dos << ", ";
    }
    // of << "\n";
    of << "],\n";

  }
  // of << "];\n";
  of << "])\n";


  //of << "x=1:" << _device_n_dofs << ";\n";
  //of << "pcolor(x, energy, log(abs(LDOS))), shading flat\n";
  
  //of << "pcolor(x, energy, abs(LDOS)), shading flat\n";
  //of << "ylabel('Energy')\n";
  //of << "xlabel('x')\n";

  of << "plt.contourf(x, energy, np.abs(LDOS), levels=10)\n";
  of << "plt.xlabel('x')\n";
  of << "plt.ylabel('Energy [eV]')\n";
  of << "plt.colorbar()\n";
  of << "plt.show()";
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
    potmodel->get_solution(elem, efermi_id, e_qfermi, elem->vertex_average());
    potmodel->get_solution(elem, hfermi_id, h_qfermi, elem->vertex_average());

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

   //In theory, size of "field" has to be consistent with _device_n_dofs, which is # of nodes
  //  unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_vars = _ext_module->get_number_of_bands();
   field.resize(_device_n_dofs * n_vars);

   //unsigned int n_vars = _ext_module->get_H_dim();
   //field.resize(n_vars); //FOR NOW THIS IS BEST SOLUTION

   switch (_which_integration)
   {
     case INTDENSITYEL:
     {

       if (get_option("quasi_equilibrium", false))
       {
         std::vector<double> Ec(_device_n_dofs*n_vars);
         std::vector<double> Ev(_device_n_dofs*n_vars);
         std::vector<double> muN(_device_n_dofs*n_vars);
         std::vector<double> muP(_device_n_dofs*n_vars);

         get_mu_and_bands(Ec, Ev, muN, muP);

         _libnegf->quasi_equilibrium_density(field, "el", Ec, Ev, muN, muP);

         error = 0.0;

         for (unsigned int i=0; i < field.size(); i++)
         {
           error += field[i];
         }

         error /= _device_n_dofs;

       }
       else
       {
         _libnegf->density(field, "el");


         error = 0.0;

         for (unsigned int i=0; i < field.size(); i++)
         {
           error += field[i];
         }

         error /= _device_n_dofs;
       }
       break;
     }

     case INTDENSITYHL:
     {
       if (get_option("quasi_equilibrium", false))
       {
         std::vector<double> Ec(_device_n_dofs*n_vars);
         std::vector<double> Ev(_device_n_dofs*n_vars);
         std::vector<double> muN(_device_n_dofs*n_vars);
         std::vector<double> muP(_device_n_dofs*n_vars);

         get_mu_and_bands(Ec, Ev, muN, muP);

         _libnegf->quasi_equilibrium_density(field, "hl", Ec, Ev, muN, muP);

         error = 0.0;

         for (unsigned int i=0; i < field.size(); i++)
         {
            error += field[i];
         }

         error /= _device_n_dofs;
       }
       else
       {

         _libnegf->density(field, "hl");

         error = 0.0;

         for (unsigned int i=0; i < field.size(); i++)
         {
           error += field[i];
         }

         error /= _device_n_dofs;
         }
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

   // WARNING: for now it works only when calling Currents
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
  if ( _hamil_type == "etb")  
  {
    transfer_density_etb(density, particle);
  }
  else if ( _hamil_type == "efa")
  {
    transfer_density_efa(density, particle);
  }  
  else 
  {
    Messages::error("undefined hamiltonian type "+_hamil_type);
  }
}


void
Negf::transfer_density_efa(const std::vector<double>& density, const std::string& particle)
{

  // compute total number of n_vars n_dofs
  // unsigned int n_vars = _sys_H->n_vars();
  unsigned int n_vars = _ext_module->get_number_of_bands();
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
        // std::cout<<band*n_dofs+n<<" dens: "<<dens<<std::endl;


        double val = equ * dens / connectivity[dof_indices_qdens[n]];

        qdens.add( dof_indices_qdens[n], val);
      }

      //for (unsigned int n = 0; n < elem->n_nodes(); n++)
      //std::cout<<connectivity[ dof_indices_qdens[n]]<<std::endl;
    }
  }

  qdens.close();
}

void
Negf::transfer_density_etb(const std::vector<double>& density, const std::string& particle)
{

  // compute total number of n_vars n_dofs
  unsigned int n_vars = _ext_module->get_number_of_bands();
  // unsigned int n_dofs = _qc_n_dofs[_quantum_contacts.size()-1];
  unsigned int N_atoms = _device_n_dofs;
  NumericVector<Number>& qdens = *_qdens_sys->solution;
  //qdens.zero();

  int particle_id;

  if (particle == "el") particle_id = _qdens_sys->variable_number("edens");
  if (particle == "hl") particle_id = _qdens_sys->variable_number("hdens");


  DofMap& dof_map_qdens = _qdens_sys->get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;

  // Calculate atomic charge by summing the density from orbitals of the same atom
  double charge;
  std::vector<double> atomic_charges(N_atoms);
  for (unsigned int i=0; i<N_atoms; i++)
  {
    charge = 0;
    for (unsigned int k=0; k<n_vars; k++)
      charge += density[i*n_vars + k];

    atomic_charges[i] = charge;
  }

  // compute connectivity of each node
  // we need the connectivity of the nodes to not double count
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

  // the cutoff distance for neighbouring atoms in A
  double cutoff = _ext_module->get_options().get_option("projection_length", 5.0);

  MeshBase::const_element_iterator el = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
  
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    dof_map_qdens.dof_indices(elem, dof_indices_qdens, particle_id);

    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      //if (qdens(dof_indices_e[n]) < 0)
      {
        auto density = project_density(elem, elem->point(n), atomic_charges, cutoff);
        density = density / connectivity[dof_indices_qdens[n]];
        qdens.add(dof_indices_qdens[n], density);
      }
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

  opt.Estep_coarse = sol_opt.get_option("Estep_coarse",1.0);

  opt.Np_n.resize(2, 0);

  sol_opt.get_option("Np_contour", opt.Np_n);

  opt.delta = sol_opt.get_option("delta",2e-3);

  opt.deltaE =  sol_opt.get_option("deltaE",opt.delta/2.0);

  int Np = 2 * opt.n_kT/opt.deltaE;

  opt.Np_real = sol_opt.get_option("Np_real",Np);

  opt.DEc = sol_opt.get_option("deltaEc", 0.25);

  opt.DEv = sol_opt.get_option("deltaEv", 0.25);

  opt.n_poles = sol_opt.get_option("Npoles", 0);

  opt.verbosity = get_option("verbosity", this->verbose());
  opt.verbosity = sol_opt.get_option("verbosity", opt.verbosity);
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
    std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));
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
    std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));
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
          unique_ptr<const Elem> elside = elem->build_side_ptr(side); 		
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

  // setup initial permutation vector as identity
  // the vector runs only on the device region where the
  // dof reordering is performed
  _perm.resize(_device_n_dofs, 0);
  for (unsigned int i = 0; i < _device_n_dofs; i++)
    _perm[i] = i;

  // Number of dofs of the first contact
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


  if (nPL < _quantum_contacts.size())
  {
    nPL = _quantum_contacts.size();
  }

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

  auto bl_it = std::unique(_end_blocks.begin(), _end_blocks.end());
  _end_blocks.resize(distance(_end_blocks.begin(), bl_it));
  nPL = _end_blocks.size();


  ostringstream os;
  os << "# DOFs       : " << n_dofs_total << "\n";
  os << "# PLs        : " << nPL << "\n";
  os << "mean PL size : " << n_dofs_total/nPL << "\n";
  Messages::info(os.str());


  if (verbose() > 3)
  {
    cerr << "PLs (block end indices):\n";
    for (auto&& b : _end_blocks)
      cerr << b << " ";
      cerr << endl;
  }



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

  /*
  for (auto&& p : _perm)
    cerr << p+1 << " ";
  cerr << endl;

  for (auto&& p : _inv_perm)
    cerr << p+1 << " ";
  cerr << endl;
  */
}



//Approach to Laplace problem
void
Negf::reorder_assemble(void)
{
  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  libMesh::DofMap& dof_map = _sys->get_dof_map();

  unsigned int n_vars = _sys->n_vars();

  libMesh::FEType fe_type = dof_map.variable_type(0);

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

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
        if (_device_n_dofs < dof_indices[n]+1)
          _device_n_dofs = dof_indices[n]+1;

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
    auto       el     = mesh.active_subdomain_set_elements_begin(qc_ids);
    const auto end_el = mesh.active_subdomain_set_elements_end(qc_ids);

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

  std::unique_ptr<FEBase> fe( FEBase::build(dim, FEType() ));

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

  for ( ; it != it_end ; ++it) 
  {
    const Elem* elem = *it;

    if (elem->subdomain_id() == qc_id)
    {
      fe->reinit(elem);

      for (unsigned int qp=0; qp<q_point.size(); qp++)
      {
        std::pair<const Elem*, Point> pair = qc->project_on_boundary(elem, q_point[qp]);

        av_V += pot_model.get_potential(pair.first, pair.second) * elem->volume()/volume;
        av_mu_n += mue_model.get_potential(pair.first, pair.second) * elem->volume()/volume;
        av_mu_p += muh_model.get_potential(pair.first, pair.second) * elem->volume()/volume;

      }
    }
  }
}

void
Negf::get_mu_and_bands(std::vector<double>& Ec, std::vector<double>& Ev,
                       std::vector<double>& muN, std::vector<double>& muP)
{
  SimulationInterface* model;

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

  //Inefficient, it performs 2x4 loops on the elements instead of 2x1. But code more clear
  Ec = get_ordered_solution(model, "Ec");
  Ev = get_ordered_solution(model, "Ev");
  muN = get_ordered_solution(model, "eQFermi");
  muP = get_ordered_solution(model, "hQFermi");

  return;
}

std::vector<double>
Negf::get_ordered_solution(SimulationInterface* model, const std::string& var)
{
  ID ID = model->get_solution_id(var);
  unsigned int n_vars = _sys_H->n_vars();

  std::vector<double> solution;
  solution.reserve(_device_n_dofs*n_vars);

  std::vector<double> tot_vals;
  tot_vals.resize(_device_n_dofs, 0.0);

  libMesh::DofMap& dof_map = _sys->get_dof_map();
  std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  assert(el != end_el);


  for (; el != end_el ; ++el )
  {
    const Elem* elem = *el;
    std::vector<double> values(elem->n_nodes());
    std::vector<Point> p(elem->n_nodes());
    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
         p[i] = elem->point(i);

    model->get_solution(elem, ID, values, p);

    dof_map.dof_indices(elem, dof_indices);

    for (size_t i = 0; i < elem->n_nodes(); ++i)
      tot_vals[dof_indices[i]] = values[i];
  }


  for (size_t i = 0; i < _device_n_dofs; ++i)
  {
    size_t index = _perm[i];

    for (size_t j = 0; j < n_vars; ++j)
      solution.push_back(tot_vals[index]);
  }


  return solution;
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
    model = nullptr;
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


double
Negf::project_density(const Elem* elem, const Point& point, const std::vector<double>& atomic_charges, double cutoff)
{
  double density(0);

  AtomisticStructure* as = _ext_module->get_atomistic_structure();

  const double scale = _ext_module->get_atomistic_structure()->get_scale();
  const double sigma = cutoff;
  const double sigma2 = 2.0*sigma*sigma;

  int N_atoms = atomic_charges.size();

  // the point in Angstrom
  Point coord(point);
  coord *= get_mesh_units() / 1e-10;

  // normalization/scale factor
  double normalization = 1.0;

  unsigned int dim = get_mesh().mesh_dimension();


  switch (dim)
  {
//      case 1:
//        normalization = 1e8;
//        break;

//      case 2:
//        normalization = 1e16 / (2.0 * M_PI * sigma);
//        break;

    default:
    {
      double tmp = 1 / (2.0 * M_PI * sigma * sigma);
      normalization = 1e24 * sqrt(tmp * tmp * tmp);
      break;
    }
  }

  int index = as->find_nearest_atom(elem, point, 2 * cutoff / 10.0);
  AtomisticBasis::neighbor_iterator it(
       as->neighbors_begin(index, cutoff));
  const AtomisticBasis::neighbor_iterator end(
       as->neighbors_end(index, cutoff));

  for (; it != end; ++it)
  {
    const Atom* atom = *it;
    unsigned int atom_id = _inv_perm[it.atom_index()];

    // atom_id may belong to a contact, skip loop in this case
    if (atom_id > N_atoms-1) continue;

    //cerr << " " << atom_id << endl;
    Point atom_pos(atom->get_position() + it.atom_translation());

  Point delta_r = coord - atom_pos;

  switch (dim)
  {
    default:
    {
      double factor = normalization * std::exp(-delta_r.norm_sq() / sigma2);
      density += atomic_charges[atom_id] * factor;

      break;
    }
  }
    
  }

  return density;
}


void
Negf::set_kpoints(std::string solution)
{
  if (solution == "Current")
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_current"));
    ModelOptions kopts;
    kopts = it->second;
    if (kopts.get_option("reduced_BZ", true))
    {
      set_kpoints_reduced_BZ(_k_int_current);
    }
    else
    {
      set_kpoints_full_BZ(_k_int_current);
    }
  }                                                       
  if (solution == "elDensity" || solution == "hlDensity")
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration_density"));
    ModelOptions kopts;
    kopts = it->second;
    if (kopts.get_option("reduced_BZ", true))
    {
      set_kpoints_reduced_BZ(_k_int_density);
    }
    else
    {
      set_kpoints_full_BZ(_k_int_density);
    }
  }
}


void
Negf::set_kpoints_full_BZ(KspaceIntegration* k_int)
{
  bool uneven_distributed;
  Kspace* negf_kspace = k_int->get_k_space();

  //Get kpoints and weights and distribute them
  map<Point, double> kpoints = k_int->get_kpoints_and_weights();
  map<Point, double> extended_kpoints;

  map<Point, double>::iterator kb(kpoints.begin());
  const map<Point, double>::iterator kend(kpoints.end());

  double degeneracy_fact = negf_kspace->get_degeneracy_factor();
  
  for (; kb != kend; kb++)
  {
    Point kpt = kb->first;
    double w = kb->second;

    std::vector<Point> eq_points;
    negf_kspace->equivalent_points(kpt, eq_points, true);

    for (auto&& eq_pt : eq_points)  extended_kpoints[eq_pt] = w / degeneracy_fact;
  }

  map<DofField, double> global_kpoints_map = k_int->broadcast_kpoints(_k_comm, extended_kpoints);
  _local_k_indices = k_int->distribute_kpoints(_k_comm, global_kpoints_map, uneven_distributed);

  if (uneven_distributed)
  {
    std::ostringstream os;
    os << "The k-points are unevenly distributed along processes. Please change the number of processes,";
    os << " the number of k-points, or the `parallel_k_groups` option to distribute them equally.";
    throw RuntimeException(os.str());
  }

  //Extract kpoints and weights separately from the map
  std::vector<double> kweights(global_kpoints_map.size());
  _global_abs_kpoints.resize(global_kpoints_map.size());

  map<DofField, double>::iterator kp_it(global_kpoints_map.begin());
  const map<DofField, double>::iterator kp_end(global_kpoints_map.end());
  int i=0;
  for (; kp_it != kp_end; kp_it++)
  {
    DofField kp = kp_it->first;
    double w = kp_it->second;
    kweights[i] = w;
    _global_abs_kpoints[i] = kp;
    i++;
  }
    
  // libNEGF wants kpoints in fractional coordinates
  std::vector<DofField> global_frac_kpoints = transform_to_fractional_coordinates(negf_kspace, _global_abs_kpoints);

  std::vector<DofField>  equivalent_points(1, std::vector<double>(3, 0.0));
  std::vector<int> n_equivalent(1, 1);

  _libnegf->set_kpoints(global_frac_kpoints, kweights, _local_k_indices, equivalent_points, n_equivalent, 0);
}


void
Negf::set_kpoints_reduced_BZ(KspaceIntegration* k_int)
{
  bool uneven_distributed;
  Kspace* negf_kspace = k_int->get_k_space();

  //Get kpoints and weights and distribute them
  map<Point, double> kpoints = k_int->get_kpoints_and_weights();
  map<DofField, double> global_kpoints_map = k_int->broadcast_kpoints(_k_comm, kpoints);
  _local_k_indices = k_int->distribute_kpoints(_k_comm, global_kpoints_map, uneven_distributed);

  if (uneven_distributed)
  {
    std::ostringstream os;
    os << "The k-points are unevenly distributed along processes. Please change the number of processes,";
    os << " the number of k-points, or the `parallel_k_groups` option to distribute them equally.";
    throw RuntimeException(os.str());
  }

  //Extract kpoints and weights separately from the map
  std::vector<double> kweights(global_kpoints_map.size());
  _global_abs_kpoints.resize(global_kpoints_map.size());

  double degeneracy_fact = negf_kspace->get_degeneracy_factor();

  map<DofField, double>::iterator kp_it(global_kpoints_map.begin());
  const map<DofField, double>::iterator kp_end(global_kpoints_map.end());
  int i=0;
  for (; kp_it != kp_end; kp_it++)
  {
    DofField kp = kp_it->first;
    double w = kp_it->second;
    kweights[i] = w / degeneracy_fact;
    _global_abs_kpoints[i] = kp;
    i++;
  }
  
  //Find equivalent points
  std::vector<DofField> equivalent_points;
  std::vector<int> n_equivalent;
  get_equivalent_points(negf_kspace, _global_abs_kpoints, equivalent_points, n_equivalent);
  
  // libNEGF wants kpoints in fractional coordinates
  std::vector<DofField> global_frac_kpoints = transform_to_fractional_coordinates(negf_kspace, _global_abs_kpoints);
  std::vector<DofField> frac_eqv_points = transform_to_fractional_coordinates(negf_kspace, equivalent_points);

  _libnegf->set_kpoints(global_frac_kpoints, kweights, _local_k_indices, frac_eqv_points, n_equivalent, 1);
}


std::vector<DofField> 
Negf::transform_to_fractional_coordinates(Kspace* kspace, const std::vector<DofField>& abs_kpoints)
{
  std::vector<DofField> frac_kpoints(abs_kpoints.size(), DofField(3, 0.0));
  int k = 0;

  for (auto&& kp_std : abs_kpoints)
  {
    Point kp;
    for (int i=0; i<3; i++) kp(i) = kp_std[i];
    kspace->inverse_transform(kp);
    for (int i=0; i<3; i++) frac_kpoints[k][i] = kp(i);
    k++;
  }
  return frac_kpoints;
}


void
Negf::get_equivalent_points(Kspace* kspace, const vector<DofField> kpoints, 
      vector<DofField>& equiv_points, vector<int>& n_equiv)
{

  n_equiv.resize(kpoints.size());
  
  for(int i=0; i<kpoints.size(); i++)
  { 
    
    DofField kp_std = kpoints[i];
    const Point kp(kp_std[0], kp_std[1], kp_std[2]);

    // Get the equivalent points for a given kpoint
    vector<Point> points;
    kspace->equivalent_points(kp, points, true);

    //Remove the last equivalent point, which coincides with the original kpoint
    points.pop_back();

    //If there are equiv. points fill the std::vector
    if (points.size() != 0)
    {
      for (auto&& p : points)
      {
        DofField std_p(3);
        for (int j=0; j<3; j++) std_p[j] = p(j);
        equiv_points.push_back(std_p);
      }
    }

    n_equiv[i] = points.size();
  }
  return;
}


std::vector<DofField>
Negf::get_lattice_vectors(std::vector<double>& r1, std::vector<double>& r2)
{
  Tensor1 vec1_real;
  Tensor1 vec2_real;

  for (short i = 0; i < 3; i++)
    vec1_real(i + 1) = r1[i]; // /(Constants::bohr_radius / mesh_units);

  for (short i = 0; i < 3; i++)
    vec2_real(i + 1) = r2[i]; // /(Constants::bohr_radius / mesh_units);

  Tensor1 vec3_real = vectorProduct(vec1_real, vec2_real);
  vec3_real = vec3_real/norm(vec3_real);

  std::vector<DofField> lattice_vectors(3);
  // libNEGF wants them in column order
  for (int j = 0; j < 3; j++)
  {
    lattice_vectors[j].resize(3);
    lattice_vectors[j][0] = vec1_real(j+1);
    lattice_vectors[j][1] = vec2_real(j+1);
    lattice_vectors[j][2] = vec3_real(j+1);
  }

  return lattice_vectors;
}


std::vector<DofField>
Negf::get_lattice_vectors(std::vector<double>& r1)
{
  Tensor1 vec1_real;

  for (short i = 0; i < 3; i++)
    vec1_real(i + 1) = r1[i];

  Tensor1 vec2_real;
  vec2_real(1) = 0;
  vec2_real(2) = -vec1_real(3);
  vec2_real(3) =  vec1_real(2);
  
  Tensor1 vec3_real = vectorProduct(vec1_real, vec2_real);
  vec3_real = vec3_real/norm(vec3_real);

  std::vector<DofField> lattice_vectors(3);
  // libNEGF wants them in column order
  for (int j = 0; j < 3; j++)
  {
    lattice_vectors[j].resize(3);
    lattice_vectors[j][0] = vec1_real(j+1);
    lattice_vectors[j][1] = vec2_real(j+1);
    lattice_vectors[j][2] = vec3_real(j+1);
  }

  return lattice_vectors;
}


void
Negf::init_basis(void)
{
  std::vector<DofField> coordinates(3);
  for (int i = 0; i < 3; i++) coordinates[i].resize(_device_n_dofs);
  get_coordinates(coordinates);

  // Transport direction for TC is x but default in NEGF is z, so we pass it as an integer.
  // Remember: Fortran is 1-indexed
  int transport_direction = 1;
  
  // Create map from atoms (or dofs) to hamiltonian matrix
  std::vector<int> matrix_indices(_device_n_dofs);
  matrix_indices[0] = 1;
  int n_bands = _ext_module->get_number_of_bands(); // In case of ETB, this is number of orbitals
  for (int i = 1; i < _device_n_dofs; i++)
  {
    matrix_indices[i] = matrix_indices[i-1] + n_bands;
  }
  _libnegf->init_basis(coordinates, matrix_indices, _lattice_vectors, transport_direction);

  return;
}


void
Negf::get_coordinates(vector<DofField>& coordinates)
{
  if (_hamil_type == "efa")
  {
    unsigned int bands = _ext_module->get_number_of_bands();
    DofMap& dof_map = _sys_H->get_dof_map();
    std::vector<unsigned int> dof_indices;    

    MeshBase::const_element_iterator       nd     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator nd_end = this->active_local_elements_end();
    for ( ; nd != nd_end; ++nd)
    {
      const Elem* elem = *nd;
      dof_map.dof_indices(elem, dof_indices, 0);

      for (unsigned int n = 0; n < elem->n_nodes(); ++n)
      {
        // coords in nanometer (consistent with lattice vectors)
        double equ =  get_mesh_units() * 1e9;
        unsigned int id = _inv_perm[dof_indices[n]];
        // coordinates[id / bands].resize(3);
        //libnegf                    tiberCAD
        coordinates[0][id / bands] = elem->point(n)(0) * equ;
        coordinates[1][id / bands] = elem->point(n)(1) * equ;
        coordinates[2][id / bands] = elem->point(n)(2) * equ;
      }
    }
  }

  if (_hamil_type == "etb")
  {
    AtomisticStructure* as = _ext_module->get_atomistic_structure();
    unsigned int N_atoms = _device_n_dofs;
    std::vector<Atom>& atoms = as->get_structure_atoms();

    for (int i=0; i<N_atoms; i++)
    {
      // coords from Angstrom to nanometer (consistent with lattice vectors)
      double equ = 0.1;
      //libnegf                    tiberCAD
      coordinates[0][i] = atoms[i].get_position(0) * equ;
      coordinates[1][i] = atoms[i].get_position(1) * equ;
      coordinates[2][i] = atoms[i].get_position(2) * equ;

    }
  }
  return;
}


// Overrides the method in SimulationInterface.C
void
Negf::setup_mpi_comm(void)
{

  // They need to be parsed here because we need the _scattering flag already at this point
  parse_scattering_options();

  if (_scattering)
  {
    // ModelOptions& device_opts = _device->get_options();
    // if (device_opts.has_submodel("Parallel"))
    // {
    //   ModelOptions& mpi_opts = device_opts.submodels_begin("Parallel")->second;
    //   unsigned int nGroups = mpi_opts.get_option("negf_k_groups", 1);
    // }
    // else
    // {
    //   Messages::error("(negf) The inelastic NEGF code needs MPI to run. "
    //   "Please include the `Parallel` block in the `Device section`");
    // }
    unsigned int nGroups = get_solver_options().get_option("parallel_k_groups", 1);
    MPI_Comm bare_kcomm;
    MPI_Comm bare_cartcomm;

    // Initialize cartesian grid. Return cart- and k- communicators
    _libnegf->mpi_cart_init(this->get_communicator().get(), nGroups, bare_cartcomm, bare_kcomm);

    // // Initialiaze libMesh Communicators
    _k_comm.duplicate(bare_kcomm);
    _cart_comm.duplicate(bare_cartcomm);

    this->set_communicator(_cart_comm);
   
  }
  else
  {
    this->set_solver_communicator(this->get_mesh().comm());
    KspaceIntegration::create_communicator(this->get_communicator(),
                                          this->get_solver_communicator(),
                                          _k_comm);
  }
  return;
}


void
Negf::parse_scattering_options(void)
{
  ModelOptions::submodel_iterator phys_opts(get_options().submodels_begin("Physics"));
  ModelOptions& phys_block = phys_opts->second;

  _deltaz = phys_block.get_option("delta_z", 0.01);
  _cell_area = phys_block.get_option("cell_area", 1.0);

  if(phys_block.has_submodel("Scattering"))
  {
    ModelOptions::submodel_iterator scatt_it(phys_block.submodels_begin("Scattering"));
    ModelOptions& scattering_block = scatt_it->second;

    if (scattering_block.has_submodel("Elastic"))
    {
      ModelOptions::submodel_iterator el_it(scattering_block.submodels_begin("Elastic"));
      const ModelOptions::submodel_iterator el_it_end(scattering_block.submodels_end("Elastic"));
      
      for (; el_it != el_it_end; el_it++)
      {
        Interaction inter;
        const ModelOptions& el_options = el_it->second;

        std::string input_model = el_options.get_option("model", "DUMMY");
        inter.model = get_scattering_model(input_model);
        if (inter.model == DUMMY) throw InitFailedException("Elastic scattering block was added without a valid model.");

        inter.coupling = el_options.get_option("coupling", inter.coupling);
        inter.scba_niter = el_options.get_option("max_scba_iterations", inter.scba_niter);
        inter.scba_tol = el_options.get_option("scba_tolerance", inter.scba_tol);
        //orbsperatom to be filled after hamiltonian has been set up

        _interactions.push_back(inter);
      }
    }

    if (scattering_block.has_submodel("Inelastic"))
    {
      ModelOptions::submodel_iterator inel_it(scattering_block.submodels_begin("Inelastic"));
      const ModelOptions::submodel_iterator inel_it_end(scattering_block.submodels_end("Inelastic"));
      for (; inel_it != inel_it_end; inel_it++)
      {
        Interaction inter;
        const ModelOptions& inel_options = inel_it->second;

        std::string input_model = inel_options.get_option("model", "DUMMY");
        inter.model = get_scattering_model(input_model);
        if (inter.model == DUMMY) throw InitFailedException("Inelastic scattering block was added without a valid model.");

        inter.coupling = inel_options.get_option("coupling", inter.coupling);
        inter.scba_niter = inel_options.get_option("max_scba_iterations", inter.scba_niter);
        inter.scba_tol = inel_options.get_option("scba_tolerance", inter.scba_tol);

        inter.wq = inel_options.get_option("phonon_frequency", inter.wq);
        inter.eps_inf = inel_options.get_option("eps_infinity", inter.eps_inf);
        inter.eps_r = inel_options.get_option("eps_0", inter.eps_r);
        inter.D0 = inel_options.get_option("deformation_potential", inter.D0);
        inter.q0 = inel_options.get_option("screening_length", inter.q0);
        inter.tTridiagonal = inel_options.get_option("tridiagonal", inter.tTridiagonal);

        _interactions.push_back(inter);
      }
    }
    _scattering = true;
  }
  else
  {
    _scattering = false;
  }

  return;
}


void
Negf::print_interactions(void)
{
  for (int i = 0; i < _interactions.size(); i++)
  {
    Interaction inter = _interactions[i];
    std::cout << "Interaction " << i << ":"<< std::endl;
    std::cout << "Model: " << inter.model << std::endl;
    std::cout << "Coupling: " << inter.coupling << std::endl;
    std::cout << "max_scba_iterations: " << inter.scba_niter <<std::endl;
    std::cout << "scba_tolerance: " << inter.scba_tol <<std::endl;
    std::cout << "orbitals per atom: "; for (auto&& v: inter.orbsperatm) std::cout << v << " "; std::cout << std::endl;
    std::cout << "frequency: " << inter.wq <<std::endl;
    std::cout << "eps_r: " << inter.eps_r <<std::endl;
    std::cout << "eps_inf: " << inter.eps_inf <<std::endl;
    std::cout << "q0: " << inter.q0 <<std::endl;
    std::cout << "D0: " << inter.D0 <<std::endl;
    std::cout << "tridiagonal: " << inter.tTridiagonal <<std::endl;
    std::cout << std::endl;
  }
  return;
}


int
Negf::get_scattering_model(std::string input_model)
{
  int model;

  if (input_model == "dephasing_diagonal") model = DEPHDIAGONAL;
  else if (input_model == "dephasing_block") model = DEPHATOMBLOCK;
  else if (input_model == "dephasing_overlap") model = DEPHOVERLAP;
  else if (input_model == "polar_optical_phonon") model = POLAROPTICAL;
  else if (input_model == "non_polar_optical_phonon") model = NONPOLAROPTICAL;
  else if (input_model == "acoustic_phonon") model = ACOUSTICINEL;
  else model = DUMMY;
  
  return model;
}


void
Negf::setup_interactions(void)
{
  double kbT = SimulationOptions::temperature * Constants::kb;
  double cell_area = _cell_area;
  double deltaz = _deltaz;
  double elastic_tol = 10.0;
  double inelastic_tol = 10.0;

  ostringstream os;

  for (int i = 0; i < _interactions.size(); i++)
  {
    Interaction inter = _interactions[i];

    int Natoms = _device_n_dofs;
    int Norbs = _ext_module->get_number_of_bands();
    std::vector<double> coup(Natoms*Norbs, inter.coupling);

    switch (inter.model)
    {
      case DEPHDIAGONAL:
        os << "Setting local fully diagonal elastic dephasing model" << std::endl;
        Messages::info(os.str()); os.str("");
        _libnegf->set_elph_dephasing(coup, inter.scba_niter);
        if ( inter.scba_tol < elastic_tol ) elastic_tol = inter.scba_tol;
        break;
      case DEPHATOMBLOCK:
        os << "Setting local block diagonal (BD) elastic dephasing model" << std::endl;
        Messages::info(os.str()); os.str("");
        _libnegf->set_elph_block_dephasing(coup, inter.orbsperatm, inter.scba_niter);
        if ( inter.scba_tol < elastic_tol ) elastic_tol = inter.scba_tol;
        break;
      case DEPHOVERLAP:
        os << "Setting overlap mask (OM) block diagonal elastic dephasing model" << std::endl;
        Messages::info(os.str()); os.str("");
        _libnegf->set_elph_s_dephasing(coup, inter.orbsperatm, inter.scba_niter);
        if ( inter.scba_tol < elastic_tol ) elastic_tol = inter.scba_tol;
        break;
      case POLAROPTICAL:
        os << "Setting polar-optical inelastic scattering model" << std::endl;
        Messages::info(os.str()); os.str("");
        _libnegf->set_elph_polaroptical(coup, inter.wq, kbT, deltaz, inter.eps_r, 
                  inter.eps_inf, inter.q0, cell_area, inter.scba_niter, inter.tTridiagonal);
        if ( inter.scba_tol < inelastic_tol ) inelastic_tol = inter.scba_tol;
        break;
      case NONPOLAROPTICAL:
        os << "Setting non polar-optical inelastic scattering model" << std::endl;
        Messages::info(os.str()); os.str("");
        _libnegf->set_elph_nonpolaroptical(coup, inter.wq, kbT, deltaz, inter.D0, 
                  cell_area, inter.scba_niter, inter.tTridiagonal);
        if ( inter.scba_tol < inelastic_tol ) inelastic_tol = inter.scba_tol;
        break;
      
      default:
        os << "Electron-phonon model in Scattering block " << i << " is not yet supported" << std::endl;
        throw std::runtime_error(os.str());
        break;
    }

  }
  _libnegf->set_scba_tolerances(elastic_tol, inelastic_tol);
}

void
Negf::set_hamiltonians(void)
{
  _n_Hk = _local_k_indices.size();
  for (_iK=0; _iK < _n_Hk; _iK++)
  {
    int k_index = _local_k_indices[_iK];
    DofField kpoint = _global_abs_kpoints[k_index];
    for(short i=0;i<3;i++) _k_vec(i) = kpoint[i];
    setup_hamil();
    finalize();
  }
}