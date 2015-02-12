#include "AtomisticGenerator2D.h"
#include "AtomisticStructure.h"
#include "mesh.h"
#include "BondMap.h"



AtomisticGenerator2D*
AtomisticGenerator2D::create(AtomisticStructure* const as)
{
  AtomisticGenerator2D* ag =  NULL;
  ag = new AtomisticGenerator2D(as);
  return ag;
}


AtomisticGenerator2D::AtomisticGenerator2D(AtomisticStructure* const as)
{
  _dim = 2;
  _as = as;
}


AtomisticGenerator2D::~AtomisticGenerator2D(void) {};


void
AtomisticGenerator2D::build()
{
  double min_x, min_y, max_x, max_y;
  Elem* elem = NULL;
  Node* nd = NULL;

  //Common building operations
  make_conv_cell();
  make_conv_basis();

  //Check edges of segment for building structure

  elem = (*_structure_elements.begin());
  nd = elem->get_node(1);

  min_x = (*nd)(0);
  min_y = (*nd)(1);
  max_x = (*nd)(0);
  max_y = (*nd)(1);


  for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){

    elem = *it;

    for (unsigned int i = 0; i < elem->n_nodes(); i++){

      nd = elem->get_node(i);
      if ( (*nd)(0) < min_x ) min_x = (*nd)(0);
      if ( (*nd)(0) > max_x ) max_x = (*nd)(0);
      if ( (*nd)(1) < min_y ) min_y = (*nd)(1);
      if ( (*nd)(1) > max_y ) max_y = (*nd)(1);
     }
  }

  // the small shift of 1e-3 prevents from subsequently cutting away
  // atoms due to roundoff errors
  _local_origin(1) += min_x * scale + 1e-3;
  _local_origin(2) += min_y * scale + 1e-3;
  //_local_origin(3) += 0.0;

  double l1 = (fabs(max_x - min_x)) * scale;
  double l2 = (fabs(max_y - min_y)) * scale;

  //Minimum periodic direction is considered along z axis, but eventually other lenghts can be
  //specified by user in input (conventional cells along these direction are assured also in this case!!)
  //double l2 = _as->get_options().get_option("y_length", 0.0);
  double l3 = _as->get_options().get_option("z_length", 0.0);

  make_supercell( l1, l2, l3);

  //print_basis(_super_basis, "supercell.xyz");
 }



