#include "AtomisticGenerator3D.h"
#include "AtomisticStructure.h"
#include "Messages.h"

AtomisticGenerator3D*
AtomisticGenerator3D::create(AtomisticStructure* const as)
{
  AtomisticGenerator3D* ag =  NULL;
  ag = new AtomisticGenerator3D(as);
  return ag;
}


AtomisticGenerator3D::AtomisticGenerator3D(AtomisticStructure* const as)
{
    _dim = 3;
    _as = as;
}


AtomisticGenerator3D::~AtomisticGenerator3D(void) {};


void
AtomisticGenerator3D::build()
{
  double min_x, min_y, max_x, max_y, min_z, max_z;
  Elem* elem = NULL;
  Node* nd = NULL;

  Messages::debug("Calling AtomisticGenerator3D::build()");

  //<<<<<<< .mine
  //3D structure is considered cluster as far as lattice vectors 
  //are not explicitely passed
  //set_periodicity(false,false,false);

  //Common building operations
  make_conv_cell();
  make_conv_lattice();
  move_origin();
  make_conv_basis();

  //Check edges of segment for building structure

  elem = (*_structure_elements.begin());
  nd = elem->get_node(1);

  min_x = (*nd)(0);
  min_y = (*nd)(1);
  min_z = (*nd)(2);
  max_x = (*nd)(0);
  max_y = (*nd)(1);
  max_z = (*nd)(2);

  for (auto it = _structure_elements.begin(); it != _structure_elements.end(); ++it){

    elem = *it;

    for (unsigned int i = 0; i < elem->n_nodes(); i++){

      nd = elem->get_node(i);
      if ( (*nd)(0) < min_x ) min_x = (*nd)(0);
      if ( (*nd)(0) > max_x ) max_x = (*nd)(0);
      if ( (*nd)(1) < min_y ) min_y = (*nd)(1);
      if ( (*nd)(1) > max_y ) max_y = (*nd)(1);
      if ( (*nd)(2) < min_z ) min_z = (*nd)(2);
      if ( (*nd)(2) > max_z ) max_z = (*nd)(2);

    }
  }


  // the small shift of 1e-3 prevents from subsequently cutting away
  // atoms due to roundoff errors
  _local_origin(1) += min_x * scale + 1e-3;
  _local_origin(2) += min_y * scale + 1e-3;
  _local_origin(3) += min_z * scale + 1e-3;


  double l1 = (fabs(max_x - min_x)) * scale;
  double l2 = (fabs(max_y - min_y)) * scale;
  double l3 = (fabs(max_z - min_z)) * scale;

  make_supercell( l1, l2, l3);

  //print_basis(_super_basis, "supercell.xyz");

  //std::cout << "Period is " << _period << std::endl;

  Messages::debug("AtomisticGenerator3D::build() done ");

}


