// $Id$

#include "Embracing.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"

#include "mesh.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"


using namespace std;


Embracing::Embracing(SimulationInterface* outer,
    SimulationInterface* inner)
  : _outer(outer),
    _inner(inner),
    _lambda(0.0),
    _is_empty(true)
{
  if ((outer == NULL) || (inner == NULL))
    throw InitFailedException("Embracing needs valid SimulationInterfaces");
}


Embracing::~Embracing(void)
{
  _elem_list.clear();
}



void
Embracing::init(const ModelOptions& options)
{
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
  // assign to all elements a weight which tells the distance from the boundary,
  // check neighbours of newly added elements etc etc

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
      for (int s = 0; s < elem->n_sides(); s++)
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

      for (int s = 0; s < elem->n_sides(); s++)
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
      for (int s = 0; s < elem->n_sides(); s++)
      {
        const Elem* neighbour =  elem->neighbor(s);
        if (neighbour != NULL)
          if ((!in.contains_element(neighbour)) && out.contains_element(neighbour))
            _sides.insert(ElementSide(elem, s));
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
      elem->set_refinement_flag(Elem::DO_NOTHING);
    else
      elem->set_refinement_flag(Elem::INACTIVE);
  }
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
  SimulationEnvironment& in = _inner->get_environment();
  const Mesh& mesh = in.get_mesh();

  DataOutput data_output(mesh, "gmv");

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

  data_output.write_cell_data("embracing", results, names);

}



void
Embracing::calculate_mixing(void)
{
}
