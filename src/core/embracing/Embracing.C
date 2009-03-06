// $Id$

#include "Embracing.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"
#include "Control.h"

#include "mesh.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"


using namespace std;


unsigned int
Embracing::_counter = 0;


Embracing::Embracing(SimulationInterface* outer,
    SimulationInterface* inner)
  : _outer(outer),
    _inner(inner),
    _lambda(0.0),
    _is_empty(true),
    _do_plot(false)
{
  if ((outer == NULL) || (inner == NULL))
    throw InitFailedException("Embracing needs valid SimulationInterfaces");

  _counter++;
}


Embracing::~Embracing(void)
{
  _elem_list.clear();
}



void
Embracing::init(const ModelOptions& options)
{
  _do_plot = options.get_option("plot_embracing_regions", false);

  SimulationEnvironment& in = _inner->get_environment();
  double x0 = (in.get_device()).get_mesh_units();
  _lambda = options.get_option("embracing_length", _lambda) / x0;
  if (_lambda > 0.0)
  {
    generate_embracing_region();
    prepare_for_solve();
    calculate_mixing();
    plot();
    reactivate_all_elements();
    _is_empty = false;
  }
}


void
Embracing::generate_embracing_region(void)
{
  // first we have to find the boundary
  find_boundary();

  // idea: begin from boundary elements and look at all neighbours,
  // assign to all elements a weight which tells the distance from the
  // boundary, check neighbours of newly added elements etc etc

  SimulationEnvironment& in = _inner->get_environment();

  set<const Elem*> elems1;
  set<const Elem*> elems2;
  set<const Elem*>* new_ptr = &elems1;
  set<const Elem*>* old_ptr = &elems2;

  const map<const Elem*, double>::iterator list_end(_elem_list.end());

  // first the boundary elements, they get weight 0
  {
    set<ElementSide>::iterator it(_sides.begin());
    const set<ElementSide>::iterator end(_sides.end());

    for ( ; it != end; ++it)
      _elem_list[it->first] = 0.0;


    for (it = _sides.begin(); it != end; ++it)
    {
      const Elem* elem = it->first;
      const Elem* neighbour;
      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        if ((s != it->second) && ((neighbour = elem->neighbor(s)) != NULL))
        {
          if ((_elem_list.find(neighbour) == list_end) && 
              (in.contains_element(neighbour)))
          {
            double w = (neighbour->centroid() - elem->centroid()).size();
            if (w < _lambda)
            {
              _elem_list[neighbour] = w;
              new_ptr->insert(neighbour);
            }
          }
        }
      }
    }
  }

  while (new_ptr->size() > 0)
  {
    set<const Elem*>* tmp = old_ptr;
    old_ptr = new_ptr;
    new_ptr = tmp;
    new_ptr->clear();

    set<const Elem*>::iterator it(old_ptr->begin());
    const set<const Elem*>::iterator end(old_ptr->end());
    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      const Elem* neighbour;

      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        if (((neighbour = elem->neighbor(s)) != NULL) &&
            (_elem_list.find(neighbour) == list_end) &&
            (in.contains_element(neighbour)))
        {
          double w = (neighbour->centroid() - elem->centroid()).size();
          w += _elem_list[elem];
          if (w < _lambda)
          {
            _elem_list[neighbour] = w;
            new_ptr->insert(neighbour);
          }
        }
      }
    }
    old_ptr->clear();
  }
}



void
Embracing::find_boundary(void)
{

  SimulationEnvironment& out = _outer->get_environment();
  SimulationEnvironment& in = _inner->get_environment();

  const Mesh& mesh = in.get_mesh();
  unsigned int dim = mesh.mesh_dimension();

  // loop over the elements of the 'inner' simulation
  MeshBase::const_element_iterator el(mesh.elements_begin());
  const MeshBase::const_element_iterator end(mesh.elements_end());
  for ( ; el != end; ++el)
  {
    const Elem* elem = *el;

    if (in.contains_element(elem))
    {
      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        const Elem* neighbour =  elem->neighbor(s);
        if (neighbour != NULL)
          if ((!in.contains_element(neighbour)) &&
              out.contains_element(neighbour))
            _sides.insert(ElementSide(elem, s));
      }
    }
  }
}



bool
Embracing::is_in_embracing_region(const Elem* elem)
{
  bool ans = false;
  if (_elem_list.find(elem) != _elem_list.end())
    ans = true;

  return ans;
}
 


void
Embracing::find_inner_boundary(void)
{

  SimulationEnvironment& in = _inner->get_environment();

  const Mesh& mesh = in.get_mesh();
  unsigned int dim = mesh.mesh_dimension();


  // loop over the elements of the embracing region
  MeshBase::const_element_iterator el(mesh.active_elements_begin());
  const MeshBase::const_element_iterator end(mesh.active_elements_end());
  for ( ; el != end; ++el)
  {
    const Elem* elem = *el;

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      const Elem* neighbour =  elem->neighbor(s);
      if (neighbour != NULL)
      {
        ElementSide side(elem, s);
        const Elem* neighbour =  elem->neighbor(s);
        if (!_sides.count(side) && !_elem_list.count(neighbour))
          _inner_sides.insert(ElementSide(elem, s));
      }
    }
  }
}



void
Embracing::prepare_for_solve(void)
{
  SimulationEnvironment& in = _inner->get_environment();
  Mesh& mesh = in.get_mesh();

  const map<const Elem*, double>::iterator list_end(_elem_list.end());

  MeshBase::element_iterator el(mesh.elements_begin());
  const MeshBase::element_iterator end(mesh.elements_end());
  for ( ; el != end; ++el)
  {
    Elem* elem = *el;

    if (_elem_list.find(elem) != list_end)
    {
      elem->set_refinement_flag(Elem::DO_NOTHING);
    }
    else
      elem->set_refinement_flag(Elem::INACTIVE);
  }

  // find the inner bopundary for boundary conditions
  find_inner_boundary();

  // create a probably unique name
  ostringstream os;
  os << "Embracing" << _counter;
  EquationSystems& eq = in.get_device().get_equation_systems();
  _laplace = &eq.add_system<LaplaceEq>(os.str());
  if (_laplace == NULL)
  {
    string msg("Could not create solver for Embracing.\n"
        "Did you perhaps name a model \'");
    msg += os.str() + "\'?";
    throw InitFailedException(msg);
  }

  _laplace->_emb = this;

  // our variable
  _laplace->add_variable("u");
  _laplace->init();
}


void
Embracing::reactivate_all_elements(void)
{
  SimulationEnvironment& in = _inner->get_environment();
  Mesh& mesh = in.get_mesh();

  const map<const Elem*, double>::iterator list_end(_elem_list.end());

  MeshBase::element_iterator el(mesh.elements_begin());
  const MeshBase::element_iterator end(mesh.elements_end());
  for ( ; el != end; ++el)
  {
    Elem* elem = *el;
    elem->set_refinement_flag(Elem::DO_NOTHING);
  }
}




void
Embracing::plot(void)
{
  // we return immediately if nothing is to be printed
  if (!_do_plot) return;

  SimulationEnvironment& in = _inner->get_environment();
  const Mesh& mesh = in.get_mesh();

  DataOutput data_output(mesh, (_inner->get_control()).get_output_format());
  data_output.set_output_directory((_inner->get_control()).get_output_dir());

  vector<double> results(mesh.n_active_elem());
  vector<string> names(1, "weight");

  MeshBase::const_element_iterator it =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh.active_elements_end();

  unsigned int elem_number = 0;
  for ( ; it != end; ++it, elem_number++)
  {
    const Elem* elem = *it;
    results[elem_number] = _elem_list[elem];
  }

  vector<double> sol;
  vector<string> solname;
  _laplace->build_nodal_results(sol, solname);

  ostringstream os;
  os << "Embracing" << _counter;
  data_output.write_nodal_data(os.str(), sol, solname);
}



void
Embracing::calculate_mixing(void)
{
  _laplace->solve();
}


  

void
Embracing::assembly(LaplaceEq& system)
{
  const unsigned int var = system.variable_number("u");
  
  DofMap& dof_map = system.get_dof_map();
  
  FEType fe_type = dof_map.variable_type(var);

  const Mesh& mesh = system.get_mesh();
  unsigned int dim = mesh.mesh_dimension();

  AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, SECOND);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  std::vector<unsigned int> dof_indices;
  
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  MeshBase::const_element_iterator el(mesh.active_elements_begin());
  const MeshBase::const_element_iterator end_el(mesh.active_elements_end());

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    
    dof_map.dof_indices (elem, dof_indices);  
    const unsigned int n_dofs   = dof_indices.size();
    
    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    for (int qp = 0; qp < qrule.n_points(); qp++)
      for (int i = 0; i < n_dofs; i++)
        for (int j = 0; j < n_dofs; j++)
          Ke(i, j) += JxW[qp] * dphi[i][qp] * dphi[j][qp];
    

    // loop over the sides for boundary conditions
    // NOTE we use penalty-method here
    for (int s = 0; s < elem->n_sides(); s++)
    {
      bool bd = false;
      double val = 0.0;
      ElementSide side(elem, s);
      if (_sides.count(side))
        bd = true;
      else if (_inner_sides.count(side))
      {
        bd = true;
        val = 1.0;
      }

      if (bd)
      {
        for (unsigned int i = 0; i < elem->n_nodes(); i++)
        {
          if (elem->is_node_on_side(i, s))
          {
            //Ke.condense(i, i, 1e6 * val, Fe);
            Ke(i, i) += 1e6;
            Fe(i) += 1e6 * val;
          }
        }

      }
    }

    // apply dof constraints
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices); 
  }
}



void
Embracing::LaplaceEq::build_nodal_results(vector<double>& results,
    vector<string>& legend)
{
  legend.resize(1);
  legend[0] = "u";

  const unsigned int s = number();
  const unsigned int var = variable_number("u");
  const Mesh& mesh = get_mesh();

  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  unsigned int number_of_points = 0;
  for ( ; nd != nd_end; ++nd)  number_of_points++;

  results.resize(number_of_points, 0.0);
  
  MeshBase::const_element_iterator it(mesh.active_local_elements_begin());
  const MeshBase::const_element_iterator
    end(mesh.active_local_elements_end());

  vector<unsigned int> dof_indices;
  DofMap& dof_map = get_dof_map();

  for ( ; it != end; ++it)
  { 
    const Elem* elem = *it;

    dof_map.dof_indices(elem, dof_indices); 

    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      unsigned int id =  elem->node(n);
      results[id]  =  (*solution)(dof_indices[n]);       
    }
  }
}


double
Embracing::get_mixing_coefficient(const Elem* elem, const Point& p)
{
  double mixing = 1.0;

  if (!_is_empty)
  {
    if (_elem_list.count(elem))
    {
      mixing = 0.0;

      LaplaceEq& system = *_laplace;

      const unsigned int var = system.variable_number("u");
      DofMap& dof_map = system.get_dof_map();
      FEType fe_type = dof_map.variable_type(var);

      const Mesh& mesh = system.get_mesh();
      unsigned int dim = mesh.mesh_dimension();

      vector<unsigned int> dof_indices;
      dof_map.dof_indices(elem, dof_indices); 

      AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
      const vector<vector<double> >& phi = fe->get_phi();

      vector<Point> points(1, p);
      FEInterface::inverse_map(dim, fe_type, elem,
          vector<Point>(1, p), points);
      fe->reinit(elem, &points);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        mixing += (*(system.solution))(dof_indices[n]) * phi[n][0];

    }

    /*
    if (!_elem_list.count(elem))
    {
      // perhaps the parent is?
      const Elem* parent = elem->parent();

      while ((parent != NULL) && (!env.contains_element(parent)))
        parent = parent->parent();

      el = parent; // is NULL if no parent

    }
    */
  }

  return mixing;
}
