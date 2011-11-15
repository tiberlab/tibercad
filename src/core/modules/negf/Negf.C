// $Id: Negf.h 2964 2011-10-10 20:34:57Z maufder $

#include "Negf.h"
#include "Device.h"
#include "QuantumContact.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"
#include "Messages.h"
#include "Boundary.h"
#include "InitFailedException.h"


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

// This is needed in order to create the shared module library
// The first string is the class name of the object to be created,
// the second one is the name of the module as it should be referred
// in the input file (the Makefile defines MODULE_NAME, which can be used here).
TIBER_MODULE(Negf, MODULE_NAME)

Negf* Negf::static_this;

Negf::Negf(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
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

  std::vector<ID> bd_ids;

  for (ID k=0; k<_contacts.size(); k++)
  {
    _device->get_boundary_region_ids(_contacts[k], bd_ids);
  }

  reorder();

}
void
Negf::do_solve(void)
{

}

void
Negf::parse_options(void)
{

  if (get_options().has_submodel("Boundary"))
  {
    _contacts.clear();

    ModelOptions::submodel_iterator it(get_options().submodels_begin("Boundary"));
    ModelOptions::submodel_iterator end(get_options().submodels_end("Boundary"));
    for( ;it!=end; ++it)
    {
      const ModelOptions& opts = it->second;
      _contacts.push_back(opts.get_name());
    }

  }
  else
  {
    throw InitFailedException("ERROR: no quantum contact");
  }

}

void
Negf::reorder(void)
{
  MeshBase& mesh = get_mesh();

  create_equation_system("linear");

  _sys = &get_equation_system<TiberLinearSystem>();

  _sys->add_variable("u", FIRST, LAGRANGE);

  _sys->attach_assemble_function(reorder_assemble);

  _sys->init();

  _sys->solve();

  const unsigned int dim = mesh.mesh_dimension();
  const NumericVector<Number>& solution = _sys->get_solution_vector();
  const DofMap& dof_map = _sys->get_dof_map();
  ID u_var = _sys->variable_number("u");

  std::vector<unsigned int> dof_indices_u;

  unsigned int sol_size = solution.size();


  std::vector<unsigned int> perm(sol_size,0);

  std::cerr<<"sol size "<<sol_size<<std::endl;

  perm.clear();
  perm.resize(sol_size, 0);

  for (unsigned int i = 0; i < sol_size; i++)
    perm[i]=i;

  //std::cerr<<"SORTING "<<std::endl;

  std::sort(perm.begin(), perm.end(), compare);

  _permu.clear();
  _permu.resize(sol_size, 0);
  for (unsigned int i = 0; i < sol_size; i++)
    _permu[perm[i]]= i;

  ID assign = 0;
  std::vector<unsigned int> temp(sol_size,sol_size);
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices_u, u_var);

    const unsigned int n_dofs = dof_indices_u.size();

    for (unsigned int i = 0; i < n_dofs; i++)
    {
      if (temp[dof_indices_u[i]] == sol_size)
      {
        //std::cerr<<"dof indices: "<<dof_indices_u[i]<<"  "<<solution(dof_indices_u[i])<<"  "<<_permu[dof_indices_u[i]]<<" "<<solution(perm[dof_indices_u[i]])<<"  "<<std::endl;
        //_permu[assign] = dof_indices_u[i];
        temp[dof_indices_u[i]] = assign;
        assign++;
      }
    }

  }

}

bool
Negf::compare(ID i, ID j)
{
  return static_this->do_compare(i, j);
}

bool
Negf::do_compare(ID i, ID j)
{
  const NumericVector<Number>& solution = _sys->get_solution_vector();
  return ( solution(i) < solution(j) );
}

void
Negf::reorder_assemble(EquationSystems& es, const std::string& system_name)
{
  static_this->do_reorder_assemble(es, system_name);
}

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

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    // Messages::info("elem init");

    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);
    Ke.zero();

    // Messages::info("elem init");
    // Messages::info("Ke");

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      for (unsigned int i=0; i<phi.size(); i++)
        for (unsigned int j=0; j<phi.size(); j++)
          Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

    // Messages::info("Fe");
    Fe.zero();

    // Messages::info ("BCs");

    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      const ElementSide elside(elem->top_parent(), side);

      Boundary* bd = _env->get_boundary(elside);

      if (bd != NULL && bd->get_name()==_contacts[0])
      {

        for(unsigned int n = 0; n< n_dofs; ++n)
        {

          if (elem->is_node_on_side(n,side))
          {

            for (unsigned int nc = 0; nc < n_dofs; nc++)
              Ke(n,nc) = 0.0;

            Ke(n,n) = 1.0;

            Fe(n) = 0.0; //fix potential

          }
        }
      }

      if (bd != NULL && bd->get_name()==_contacts[1])
      {
        for(unsigned int n = 0; n< n_dofs; ++n)
        {

          if (elem->is_node_on_side(n,side))
          {
            for (unsigned int nc = 0; nc < n_dofs; nc++)
              Ke(n,nc) = 0.0;

            Ke(n,n) = 1.0;

            Fe(n) = 1.0; //fix potential

          }
        }
      }


      //Messages::info("BCs");

    }

    dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

    //Messages::info ("matrix insertion");

    _sys->matrix->add_matrix (Ke, dof_indices);
    _sys->rhs->add_vector    (Fe, dof_indices);

    //Messages::info("matrix insertion");

  }

}



void
Negf::do_setup_solution_variables(void)
{
  declare_solution(ReorderPotential, REAL, NODES, "1");
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
    ID u_var = _sys->variable_number("u");
    FEType fe_type = _sys->variable_type(u_var);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &p);

    std::vector<unsigned int> dof_indices_u;

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    const unsigned int n_dofs = dof_indices_u.size();
    unsigned int np = p.size();

    //std::cout<<"elem id "<<elem->id()<<std::endl;

    for (unsigned int n = 0; n < np; n++)
    {
      double u  = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        //std::cout<<"dof ind "<<dof_indices_u[i]<<std::endl;
        u += phi[i][n] * solution(dof_indices_u[i]); ;
      }

      values[ReorderPotential][n] = u;

    }
  }

}

void
Negf::test_project_on_boundary(void)
{
  QuantumContact* qc = _device->get_quantum_contact(_contacts[0]);

  MeshBase& mesh = get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  AutoPtr<FEBase> fe( FEBase::build(dim, FEType() ));

  QGauss qrule(dim, FIRST); // regola di ordine 0 perche in questo modo prendo la normale del centroide

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
      {//qp

        //std::cerr<<"qpnt: "<<q_point[qp](0)<<" "<<q_point[qp](1)<<" "<<q_point[qp](2)<<std::endl;

        std::pair<const Elem*, Point> pair = qc->project_on_boundary(elem, q_point[qp]);

        Point pt = pair.second;
        const Elem* sidelem = pair.first;

        //std::cerr<<"Pt: "<<pt(0)<<" "<<pt(1)<<" "<<pt(2)<<std::endl;

        //std::cerr<<"side elem id: "<<sidelem->id()<<std::endl;

      }

    }
  }
}
