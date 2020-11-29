// $Id$

#include "Embracing.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"
#include "TiberCad.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"


using namespace std;



unsigned int
Embracing::_counter = 0;


Embracing::Embracing(SimulationInterface* outer,
    SimulationInterface* inner)
  : _outer(outer),
    _inner(inner),
    _lambda(0.0),
    _cutoff(0.0),
    _is_empty(true),
    _do_plot(false),
    _need_mixing(false),
    _laplace(NULL)
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
  _do_plot = options.get_option("plot_embracing_region", _do_plot);
  _need_mixing = options.get_option("calculate_mixing", _need_mixing);

  if (!_inner->has_environment() || !_outer->has_environment())
  {
    ostringstream msg;
    msg << "Could not create Embracing between \'"
        << _inner->get_name() << "\' and \'"
        << _outer->get_name() << "\'.\n\'"
        << (_inner->has_environment() ? _outer->get_name() : _inner->get_name())
        << "\' does not have a mesh.";
    throw InitFailedException(msg.str());
  }

  double x0 = _inner->get_mesh_units();
  _lambda = options.get_option("embracing_length", _lambda) / x0;
  // cutoff is read as percentage of embracing region
  _cutoff = options.get_option("cutoff", _cutoff) * _lambda;

  if (_lambda > 0.0)
  {
    generate_embracing_region();
    if (_do_plot || _need_mixing)
    {
      calculate_mixing();
      plot();
    }
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

  elem_iterator list_end(_elem_list.end());

  // first the boundary elements, they get weight 0
  {
    set<ElementSide>::iterator it(_sides.begin());
    const set<ElementSide>::iterator end(_sides.end());

    for ( ; it != end; ++it)
      _elem_list[it->elem()] = 0.0;


    for (it = _sides.begin(); it != end; ++it)
    {
      const Elem* elem = it->elem();
      const Elem* neighbour;
      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        if ((s != it->side()) && ((neighbour = elem->neighbor(s)) != NULL))
        {
          if ((find_elem(neighbour) == list_end) &&
              (in.contains_element(neighbour)))
          {
            double w = (neighbour->centroid() - elem->centroid()).norm();
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
            (find_elem(neighbour) == list_end) &&
            (in.contains_element(neighbour)))
        {
          double w = (neighbour->centroid() - elem->centroid()).norm();
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

  // loop over the elements of the 'inner' simulation
  SimulationEnvironment::ConstElemIterator el(in.elements_begin());
  const SimulationEnvironment::ConstElemIterator end(in.elements_end());
  for ( ; el != end; ++el)
  {
    const Elem* elem = *el;

    if (in.contains_element(elem))
    {
      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        const Elem* neighbour =  elem->neighbor(s);
        if (neighbour != NULL)
          if (out.contains_element(neighbour) &&
              (!in.contains_element(neighbour) ||
               !out.contains_element(elem)))
            _sides.insert(ElementSide(elem, s));
      }
    }
  }
}



void
Embracing::find_inner_boundary(void)
{

  SimulationEnvironment& in = _inner->get_environment();

  const MeshBase& mesh = in.get_mesh();


  // loop over the elements of the embracing region
  //MeshBase::const_element_iterator el(mesh.active_elements_begin());
  //const MeshBase::const_element_iterator end(mesh.active_elements_end());
  //for ( ; el != end; ++el)
  for (auto&& el : _elem_list)
  {
    const Elem* elem = el.first;

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
Embracing::plot(void)
{
  // we return immediately if nothing is to be printed
  if (!_do_plot) return;

  const MeshBase& mesh = _inner->get_mesh();

  DataOutput data_output(mesh, TiberCad::get_output_format());
  data_output.set_output_directory(TiberCad::get_output_dir());


  vector<double> sol;
  vector<string> solname;
  ostringstream os;
  os << "Embracing" << _counter;

  if (_laplace != NULL)
  {
    _laplace->build_nodal_results(sol, solname);
    data_output.write_nodal_data(os.str(), sol, solname);
  }
  else
  {
    sol.resize(mesh.n_active_elem());
    solname.resize(1);
    solname[0] = "region_ID";

    MeshBase::const_element_iterator it =
      mesh.active_elements_begin();
    const MeshBase::const_element_iterator end =
      mesh.active_elements_end();

    unsigned int elem_number = 0;
    for ( ; it != end; ++it, elem_number++)
    {
      const Elem* elem = *it;
      sol[elem_number] = elem->subdomain_id();
    }
    data_output.write_cell_data(os.str(), sol, solname);
  }

}



void
Embracing::calculate_mixing(void)
{
  // we return immediately if mixing coeffs are not required
  if (!_need_mixing) return;

  find_inner_boundary();

  if (_inner->get_mesh().mesh_dimension() == 1)
  {
    // ordered list of boundary points
    set<Point> bd_points;

    // ordered list of inner boundary points
    set<Point> inner_points;

    for (auto&& s : _sides)
    {
      bd_points.insert(s.elem()->point(s.side()));
    }

    for (auto&& s : _inner_sides)
    {
      inner_points.insert(s.elem()->point(s.side()));
    }

    // now the embracing regions are necessarily formed by the pairs
    // bd_points(i), inner_points(i)
    _segments.reserve(bd_points.size());

    auto it = bd_points.begin();
    auto it2 = inner_points.begin();
    for ( ; it != bd_points.end(); ++it, ++it2)
    {
      _segments.push_back(Interval(*it, Point(*it2 - *it)));
    }
  }
  else if (_laplace == NULL)
  {
    // find the inner bopundary for boundary conditions
    find_inner_boundary();

    // create a probably unique name
    ostringstream os;
    os << "Embracing" << _counter;
    SimulationEnvironment& in = _inner->get_environment();
    libMesh::EquationSystems& eq = in.get_device().get_equation_systems();
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

    _laplace->solve();
  }
}




void
Embracing::assembly(LaplaceEq& system)
{
  const unsigned int var = system.variable_number("u");

  libMesh::DofMap& dof_map = system.get_dof_map();

  libMesh::FEType fe_type = dof_map.variable_type(var);

  const MeshBase& mesh = system.get_mesh();
  unsigned int dim = mesh.mesh_dimension();

  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
  libMesh::QGauss qrule(dim, libMesh::SECOND);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  std::vector<unsigned int> dof_indices;

  libMesh::DenseMatrix<Number> Ke;
  libMesh::DenseVector<Number> Fe;

  const double penalty = 1e6;

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

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
      for (unsigned int i = 0; i < n_dofs; i++)
        for (unsigned int j = 0; j < n_dofs; j++)
          Ke(i, j) += JxW[qp] * dphi[i][qp] * dphi[j][qp];


    if ((_cutoff > 0.0) && (_elem_list[elem] < _cutoff))
      for (size_t n = 0; n < elem->n_nodes(); n++)
        Ke(n, n) += penalty;

    // loop over the sides for boundary conditions
    // NOTE we use penalty-method here
    for (unsigned int s = 0; s < elem->n_sides(); s++)
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
            Ke(i, i) += penalty;
            Fe(i) += penalty * val;
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

  //const unsigned int s = number();
  //const unsigned int var = variable_number("u");
  const MeshBase& mesh = get_mesh();

  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  unsigned int number_of_points = 0;
  for ( ; nd != nd_end; ++nd)  number_of_points++;

  results.resize(number_of_points, 0.0);

  MeshBase::const_element_iterator it(mesh.active_local_elements_begin());
  const MeshBase::const_element_iterator
    end(mesh.active_local_elements_end());

  vector<unsigned int> dof_indices;
  libMesh::DofMap& dof_map = get_dof_map();

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

  assert(_need_mixing == true);

  if (!_is_empty)
  {
    if (_inner->get_mesh().mesh_dimension() == 1)
    {
      // find the interval p is in
      for (auto&& s : _segments)
      {
        const Point& p1 = s.first;
        Point d = s.second;
        double t = (p - p1) * d / d.norm_sq();
        if ((t >= 0) && (t < 1))
        {
          mixing = t;
          break;
        }
      }

    }
    else
    {

      while ((elem != NULL) && (!_elem_list.count(elem)))
        elem = elem->parent(); // perhaps the parent is in region?

      // elem is now NULL if no parent
      // if non NULL, it is for sure in the embracing region
      if (elem != NULL)
      {
        mixing = 0.0;

        LaplaceEq& system = *_laplace;

        const unsigned int var = system.variable_number("u");
        libMesh::DofMap& dof_map = system.get_dof_map();
        libMesh::FEType fe_type = dof_map.variable_type(var);

        const MeshBase& mesh = system.get_mesh();
        unsigned int dim = mesh.mesh_dimension();

        vector<unsigned int> dof_indices;
        dof_map.dof_indices(elem, dof_indices);

        libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
        const vector<vector<double> >& phi = fe->get_phi();

        vector<Point> points(1, p);
        libMesh::FEInterface::inverse_map(dim, fe_type, elem,
            vector<Point>(1, p), points);
        fe->reinit(elem, &points);

        for (unsigned int n = 0; n < elem->n_nodes(); n++)
          mixing += (*(system.solution))(dof_indices[n]) * phi[n][0];

      }
    }
  }

  return mixing;
}
