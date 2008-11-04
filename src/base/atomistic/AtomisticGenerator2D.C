#include "AtomisticGenerator2D.h"
#include "AtomisticStructure.h"
#include "mesh.h"




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
	_rotation(1,1) = 1.0; _rotation(1,2) = 0.0; _rotation(1,3) = 0.0; _rotation(2,1) = 0.0; _rotation(2,2) = 1.0;
	_rotation(2,3) = 0.0; _rotation(3,1) = 0.0; _rotation(3,2) = 0.0; _rotation(3,3) = 1.0;

	_lattice_constant[0] = 0.0; _lattice_constant[1] = 0.0; _lattice_constant[2] = 0.0;
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

  _local_origin(1) = min_x * scale; _local_origin(2) = min_y * scale; _local_origin(3) = 0.0;

	double l1 = (fabs(max_x - min_x)) * scale;
	double l2 = (fabs(max_y - min_y)) * scale;

	//Minimum periodic direction is considered along z axis, but eventually other lenghts can be
	//specified by user in input (conventional cells along these direction are assured also in this case!!)
	//double l2 = _as->get_options().get_option("y_lenght", 0.0);
	double l3 = _as->get_options().get_option("z_lenght", 0.0);

	make_supercell( l1, l2, l3);

	//print_basis(_super_basis, "supercell.xyz");
 }



void AtomisticGenerator2D::passivate(void){

  std::vector<Atom> periodic_basis, hydrogens;
  unsigned int i;
  Atom tmp;
  unsigned int ** bond_map_periodic;
  Tensor1 z_period;

 z_period(1) = _period(1,3); z_period(2) = _period(2,3); z_period(3) = _period(3,3);

  for (i = 0; i < _structure_basis.size(); i++){

    tmp = _structure_basis[i];
    tmp.set_flag(1);
    periodic_basis.push_back(tmp);

    tmp=_structure_basis[i];
     tmp.set_position(tmp.get_position() - z_period);
     //tmp.position = tmp.position + z_period;
    periodic_basis.push_back(tmp);

    tmp=_structure_basis[i];
    tmp.set_position(tmp.get_position() - z_period);
    periodic_basis.push_back(tmp);

  }

  //bond_map_periodic = bond_map_gen(&(periodic_basis));
  //bond_map_periodic = bond_map_gen(periodic_basis);
  if (_bondmapobject != NULL) delete _bondmapobject;
  passivate_cluster(periodic_basis);

  _structure_basis.clear();

  for (i = 0; i < periodic_basis.size(); i++){
    if (periodic_basis[i].get_flag() == 1) _structure_basis.push_back(periodic_basis[i]);


}

};


