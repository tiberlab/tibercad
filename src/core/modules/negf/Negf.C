// $Id: Negf.h 2964 2011-10-10 20:34:57Z maufder $

#include "Negf.h"
#include "Device.h"
#include "QuantumContact.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"
#include "Messages.h"
#include "InitFailedException.h"
#include "PotentialInterface.h"

// Basic include files needed for the mesh functionality.
#include "fe.h"
#include "fe_interface.h"
// Define generic quadrature rules.
#include "quadrature.h"
#include "quadrature_gauss.h"
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

// C++ includes
#include <fstream>
#include <set>
#include <algorithm>
#include <math.h>

using namespace Constants;

TIBER_MODULE(Negf, MODULE_NAME)

Negf* Negf::static_this;

Negf::Negf(const ModelOptions& options) :
              SimulationInterface(options)
{
  _device_n_dofs = 0;
  _qc_n_dofs = 0;
}

Negf::~Negf(void)
{
}

Negf*
Negf::create(const ModelOptions& options)
{
  return static_this = new Negf(options);
}

//! The initialization
void
Negf::do_init(void)
{

  parse_options();

  _env = &get_environment();

  _device = &(_env->get_device());

  {
    _boundaries.clear();
    ID id = -1;
    SimulationEnvironment::BoundaryIterator it = _env->boundaries_begin();
    const SimulationEnvironment::BoundaryIterator end = _env->boundaries_end();
    for( ;it!=end; ++it)
    {
      _boundaries[*it] = ++id;
    }
  }
  // Prepare QuantumContact Map
  // Note: Boundary _contact_names are the same as QC names !!
  {
    _quantum_contacts.clear();
    std::map<const Boundary*, ID>::iterator it = _boundaries.begin();
    const std::map<const Boundary*, ID>::iterator end = _boundaries.end();
    for (; it != end; ++it)
    {
      QuantumContact* qc = _device->get_quantum_contact(it->first->get_name());

      _quantum_contacts[qc->get_id()] = qc;
      //std::cerr<<"_quantum_contact "<<qc->get_id()<<" "<<it->first->get_name()<<" "<<it->second<<std::endl;
      // Quantum Contacts are activated here: so dof_map come out correctly
      qc->activate_elements();
    }
  }
}

void
Negf::setup_effectivemass_hamil()
{
  // Setup a simple effective mass Hamiltonian
  ID id = create_equation_system("linear");
  _sys_H = &get_equation_system<TiberLinearSystem>(id);
  _sys_H->add_variable("phi", FIRST, LAGRANGE);

  id = create_equation_system("linear");
  _sys_S = &get_equation_system<TiberLinearSystem>(id);
  _sys_S->add_variable("phi", FIRST, LAGRANGE);

  _sys_H->attach_assemble_function(ham_assemble);

  _sys_H->init();
  _sys_S->init();

  reorder(); // dof indices reorder

  _sys_H->assemble();
}

void
Negf::do_solve(void)
{
  setup_effectivemass_hamil();
}

void
Negf::parse_options(void)
{

  _pot_module = get_option("potential_simulation","");

  _n_blocks = get_option("number_of_blocks", 10);

}

void
Negf::do_setup_solution_variables(void)
{
  declare_solution(ReorderPotential, REAL, NODES, "1");
}

void
Negf::ham_assemble(EquationSystems& es, const std::string& system_name)
{
  static_this->do_ham_assemble(es, system_name);
}


void
Negf::do_ham_assemble(EquationSystems& es, const std::string& system_name)
{
  PotentialInterface model;

  model.set_simulation(_pot_module);

  const double newconst = 0.5 * Hartree * bohr_radius/get_mesh_units() * bohr_radius/get_mesh_units();

  TensorValue<double> invMass(0.0);

  std::vector<double> V;

  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map = _sys_H->get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

  QGauss qrule(dim, THIRD);

  fe->attach_quadrature_rule(&qrule);

  AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));

  QGauss qface(dim-1, THIRD);

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseMatrix<Number> He; //
  DenseMatrix<Number> Se;

  std::vector<unsigned int> dof_indices,new_dof_indices;

  //ACTIVATE QC
  activate_quantum_contacts();

  //ITERATION OVER ACTIVE DEVICE REGION
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    //get effective mass tensor for elem
    //negfmod = <NegfModel>(get_bulk_model(elem));
    invMass(0,0) = 1.0; invMass(1,1) = 1.0;  invMass(2,2) = 1.0;

    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    //get potential from dd model
    V.resize(qrule.n_points());
    V.assign(qrule.n_points(), 0.0);

    if( model.has_simulation() && !model.get_simulation()->is_solved() )
    {
      model.get_potential(elem, q_point, V);
    }

    He.resize(n_dofs, n_dofs);
    Se.resize(n_dofs, n_dofs);
    He.zero();
    Se.zero();

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      for (unsigned int i=0; i<phi.size(); i++)
        for (unsigned int j=0; j<phi.size(); j++)
        {
          Se(i,j) += JxW[qp]* phi[i][qp] * phi[j][qp];
          He(i,j) += JxW[qp]* newconst * dphi[i][qp] * (invMass*dphi[j][qp]);
          He(i,j) += JxW[qp] * V[qp] * phi[i][qp] * phi[j][qp];
        }

    new_dof_indices.resize(n_dofs);

    for (unsigned int i=0; i< n_dofs; i++)
      new_dof_indices[i] = _permu[dof_indices[i]];


    _sys_S->matrix->add_matrix(Se, new_dof_indices);

    _sys_H->matrix->add_matrix(He, new_dof_indices);

  }
  _sys_H->matrix->print_matlab("H.m");
  _sys_S->matrix->print_matlab("S.m");
  std::cerr<<"print"<<std::endl;
}


void
Negf::get_solution_secure(const Elem *elem, std::map< ID, std::vector< double > > &values,
    const std::vector< Point > &p)
{

  if (values.count(ReorderPotential))
  {
    const unsigned int dim = get_mesh().mesh_dimension();
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
    unsigned int np = p.size();


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
  std::cerr<<"Reorder dofs"<<std::endl;

  MeshBase& mesh = get_mesh();

  ID id = create_equation_system("linear");
  _sys = &get_equation_system<TiberLinearSystem>(id);

  unsigned int n_vars = _sys_H->n_vars();

  //std::cerr<<"n vars "<<n_vars<<std::endl;

  std::vector<ID> u_vars(n_vars,0);

  for (unsigned int n=0; n< n_vars; n++)
  {
    std::ostringstream var_str;
    var_str << "u" << n;
    std::string name = var_str.str();
    _sys->add_variable(name, FIRST, LAGRANGE);
    u_vars[n] = _sys->variable_number(name);
    //std::cerr<<"var: "<<u_vars[n]<<std::endl;
  }

  _sys->attach_assemble_function(reorder_assemble);

  _sys->init();

  _sys->solve();

  std::cerr<<"Laplace solved"<<std::endl;

  const NumericVector<Number>& solution = _sys->get_solution_vector();

  const DofMap& dof_map = _sys->get_dof_map();

  std::vector<unsigned int> dof_indices_u;

  unsigned int sol_size = solution.size();

  std::vector<unsigned int> perm(sol_size,0);

  //std::cerr<<"sol size "<<sol_size<<std::endl;

  // setup permutation vector
  perm.clear();
  perm.resize(_device_n_dofs+1, 0);
  for (unsigned int i = 0; i < _device_n_dofs+1; i++)
    perm[i]=i;

  std::cerr<<"Sorting"<<std::endl;

  std::sort(perm.begin(), perm.end(), compare);

  // Invert permutation to compute _permu
  _permu.clear();
  _permu.resize(_qc_n_dofs+1, 0);
  for (unsigned int i = 0; i < _device_n_dofs+1; i++)
    _permu[perm[i]]= i;

  // Reset dofs in QC as an identity
  for (unsigned int i = _device_n_dofs+1; i < _qc_n_dofs+1; i++)
    _permu[i]= i;

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
        //std::cerr<<"dof indices: "<<dof_indices_u[i]<<"  "<<solution(dof_indices_u[i])<<"  ";
        //std::cerr<<_permu[dof_indices_u[i]]<<" "<<solution(_permu[dof_indices_u[i]])<<"  "<<std::endl;
        temp[dof_indices_u[i]] = assign;
        assign++;
      }
    }
  }

  // Create slicing 
  double slice = 1.0/_n_blocks;
  unsigned int i_slice = 0;
  _end_blocks.resize(_n_blocks);
  _end_blocks.assign(_n_blocks, -1);

  for (unsigned int i = 0; i < _device_n_dofs+1; i++)
  {
    if (solution(perm[i]) <= slice+0.005 )
    {
      _end_blocks[i_slice]++;
    }
    else
    {
      i_slice++;
      slice = (i_slice+1) * 1.0/_n_blocks;
      _end_blocks[i_slice] = _end_blocks[i_slice-1] + 1;
    }
  }
  std::cerr<<"Slice created"<<std::endl;
/*
  for (unsigned int i = 0; i < _end_blocks.size(); i++)
  {
    std::cerr<<_end_blocks[i]<<std::endl;
  }
*/
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
  const MeshBase& mesh = get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map = _sys->get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

  QGauss qrule(dim, THIRD);

  fe->attach_quadrature_rule(&qrule);

  AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));

  QGauss qface(dim-1, THIRD);

  fe_face->attach_quadrature_rule(&qface);

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

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

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

      if (bd != NULL)
      {
        ID cc = _boundaries[bd];

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
      if (_device_n_dofs < dof_indices[n]) _device_n_dofs = dof_indices[n];
  }
  //std::cerr<<"n dofs in device: "<<_device_n_dofs<<std::endl;

  // Iterates over INACTIVE QuantumContact regions
  {
    MeshBase::const_element_iterator       el     = mesh.not_active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.not_active_elements_end();

    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      // Checks whether the element is within a QuantumContact
      ID sub_id = elem->subdomain_id();
      if (_quantum_contacts.count(sub_id))
      {
        dof_map.dof_indices(elem, dof_indices);
        const unsigned int n_dofs= dof_indices.size();

        fe->reinit(elem);
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);
        Ke.zero();
        Fe.zero();

        for(unsigned int n = 0; n< n_dofs; ++n)
        {
          Ke(n,n) = penalty;
          ID cc;

          std::string qc_name = _device->get_region_name(sub_id);
          cc = _boundaries[_env->get_boundary(qc_name)];

          Fe(n) = penalty * cc;
          // Update number of dofs in quantum contact. It start from number of dofs in device region
          if (_qc_n_dofs < dof_indices[n]) _qc_n_dofs = dof_indices[n];
        }
      }
      _sys->matrix->add_matrix(Ke, dof_indices);
      _sys->rhs->add_vector(Fe, dof_indices);
    }
    //std::cerr<<"dofs in qc: "<<_qc_n_dofs - _device_n_dofs<<std::endl;
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
Negf::test_project_on_boundary(void)
{
  QuantumContact* qc = _device->get_quantum_contact(_boundaries.begin()->first->get_name());

  MeshBase& mesh = get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  AutoPtr<FEBase> fe( FEBase::build(dim, FEType() ));

  QGauss qrule(dim, FIRST); // Order 0 rule because in this way we take centroid's normal

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Point>& q_point = fe->get_xyz();

  ID qc_id =qc->get_id();

  MeshBase::const_element_iterator it= mesh.active_elements_begin();
  const MeshBase::const_element_iterator it_end  = mesh.active_elements_end();

  for ( ; it != it_end ; ++it) //loop over k space elements
  {
    const Elem* elem = *it;

    if (elem->subdomain_id() == qc_id)
    {
      fe->reinit(elem);

      for (unsigned int qp=0; qp<q_point.size(); qp++)
      {
        std::pair<const Elem*, Point> pair = qc->project_on_boundary(elem, q_point[qp]);

        Point pt = pair.second;
        const Elem* sidelem = pair.first;
      }
    }
  }
}
