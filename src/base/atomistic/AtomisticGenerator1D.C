#include "AtomisticGenerator1D.h"
#include "AtomisticStructure.h"




AtomisticGenerator1D*
AtomisticGenerator1D::create(AtomisticStructure* const as)
{
  AtomisticGenerator1D* ag =  NULL;
  std::cout << "creating 1D atomistic structure... ";
  ag = new AtomisticGenerator1D(as);
  return ag;
}


AtomisticGenerator1D::AtomisticGenerator1D(AtomisticStructure* const as)
{
  _dim = 1;
  _as = as;
}


AtomisticGenerator1D::~AtomisticGenerator1D(void)
{
}


void
AtomisticGenerator1D::build()
{

  // <<<<<<< .mine
  //1D and 2D structures are intended to be periodical
  //set_periodicity(false,true,true);

  Elem* elem = NULL;
  Node* nd = NULL;

  //Common building operations
  make_conv_cell();
  if (_as->get_options().get_option("minimal_cell",false))
  {
    minimal_conv_cell();
  }
  make_conv_lattice();
  move_origin();
  make_conv_basis();

  //Check edges of segment for building structure
  //MeshBase::node_iterator nd = _as->get_device()->get_mesh().nodes_begin();
  //    const MeshBase::node_iterator nd_end = _as->get_device()->get_mesh().nodes_end();
  elem = (*_structure_elements.begin());
  nd = elem->get_node(1);

  double edge_min, edge_max;

  edge_min = (*nd)(0); edge_max = (*nd)(0);
  nd++;

  for (auto it = _structure_elements.begin(); it != _structure_elements.end(); it++)
  {

    elem = *it;

    for (unsigned int i = 0; i < elem->n_nodes(); i++)
    {
      nd = elem->get_node(i);
      if ( (*nd)(0) < edge_min ) edge_min = (*nd)(0);
      if ( (*nd)(0) > edge_max ) edge_max = (*nd)(0);
    }
  }

  // the small shift of 1e-3 prevents from subsequently cutting away
  // atoms due to roundoff errors
  _local_origin(1) += edge_min * scale + 1e-3;
  //_local_origin(2) += 0.0;
  //_local_origin(3) += 0.0;
  double l1 = (fabs(edge_max - edge_min)) * scale;

  //Minimum periodic direction is considered along y and z axis, but eventually other lenghts can be
  //specified by user in input (conventional cells along these direction are assured also in this case!!)
  double l2 = _as->get_options().get_option("supercell_size_y", 0.0);
  double l3 = _as->get_options().get_option("supercell_size_z", 0.0);

  make_supercell( l1, l2, l3);

}


