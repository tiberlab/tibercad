#include "AtomisticGenerator3D.h"
#include "AtomisticStructure.h"
#include "mesh.h"
#include "BondMap.h"


AtomisticGenerator3D*
AtomisticGenerator3D::create(AtomisticStructure* const as)
{
  AtomisticGenerator3D* ag =  NULL;
  std::cout << "creating 3D atomistic structure... ";
  ag = new AtomisticGenerator3D(as);
  std::cout << "done" << std::endl;
  return ag;
}


AtomisticGenerator3D::AtomisticGenerator3D(AtomisticStructure* const as)
{
  _dim = 3;
                   _as = as;
	_rotation(1,1) = 1.0; _rotation(1,2) = 0.0; _rotation(1,3) = 0.0; _rotation(2,1) = 0.0; _rotation(2,2) = 1.0;
	_rotation(2,3) = 0.0; _rotation(3,1) = 0.0; _rotation(3,2) = 0.0; _rotation(3,3) = 1.0;

	_lattice_constant[0] = 0.0; _lattice_constant[1] = 0.0; _lattice_constant[2] = 0.0;
}


AtomisticGenerator3D::~AtomisticGenerator3D(void) {};


void
AtomisticGenerator3D::build()
{
  double min_x, min_y, max_x, max_y, min_z, max_z;
  Elem* elem = NULL;
  Node* nd = NULL;

#ifdef DEBUG
  std::cerr << "Calling AtomisticGenerator3D::build() " << std::endl;
#endif

  //Common building operations
  make_conv_cell();
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


  for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){

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


  std::cerr << "min_x = " << min_x << std::endl;
  std::cerr << "max_x = " << max_x << std::endl;
  std::cerr << "min_y = " << min_y << std::endl;
  std::cerr << "max_y = " << max_y << std::endl;
  std::cerr << "min_z = " << min_z << std::endl;
  std::cerr << "max_z = " << max_z << std::endl;

  _local_origin(1) = min_x * scale; _local_origin(2) = min_y * scale; _local_origin(3) = min_z * scale;

  std::cerr << "_local_origin is " << _local_origin << std::endl;

  double l1 = (fabs(max_x - min_x)) * scale;
  double l2 = (fabs(max_y - min_y)) * scale;
  double l3 = (fabs(max_z - min_z)) * scale;

  make_supercell( l1, l2, l3);

  //print_basis(_super_basis, "supercell.xyz");

  //std::cout << "Period is " << _period << std::endl;

#ifdef DEBUG
  std::cerr << "AtomisticGenerator3D::build() done " << std::endl;
#endif

}



void AtomisticGenerator3D::passivate(void){

  unsigned int ** bond_map;

  //bond_map_periodic = bond_map_gen(&(periodic_basis));
  //bond_map = bond_map_gen(_structure_basis);
  std::cout << "_bondmapobject " << _bondmapobject << std::endl;
  if (_bondmapobject != NULL) {delete _bondmapobject; _bondmapobject = new BondMap;}
  std::cout << "_bondmapobject " << _bondmapobject << std::endl;

  passivate_cluster(_structure_basis);


};
