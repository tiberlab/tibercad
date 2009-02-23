#include "AtomisticGenerator1D.h"
#include "AtomisticStructure.h"
#include "mesh.h"
#include "BondMap.h"




AtomisticGenerator1D*
AtomisticGenerator1D::create(AtomisticStructure* const as)
{
  AtomisticGenerator1D* ag =  NULL;
  std::cout << "creating 1D atomistic structure... ";
  ag = new AtomisticGenerator1D(as);
  std::cout << "done" << std::endl;
  return ag;
}


AtomisticGenerator1D::AtomisticGenerator1D(AtomisticStructure* const as)
{
std::cout << "start constructor" << std::endl;
  _dim = 1;
                   _as = as;
	_rotation(1,1) = 1.0; _rotation(1,2) = 0.0; _rotation(1,3) = 0.0; _rotation(2,1) = 0.0; _rotation(2,2) = 1.0;
	_rotation(2,3) = 0.0; _rotation(3,1) = 0.0; _rotation(3,2) = 0.0; _rotation(3,3) = 1.0;

	_lattice_constant[0] = 0.0; _lattice_constant[1] = 0.0; _lattice_constant[2] = 0.0;
std::cout << "end constructor" << std::endl;
}


AtomisticGenerator1D::~AtomisticGenerator1D(void){};


 void
 AtomisticGenerator1D::build()
 {

   //Only 1D and 2D structures are intended to be periodical
   _as->_atomistic_structure_options.is_periodical = true;

	 Elem* elem = NULL;
	 Node* nd = NULL;

//Common building operations
	 make_conv_cell();
	 make_conv_basis();

	 //unsigned int dimension = _as->get_device()->get_mesh().mesh_dimension();

	 //Check edges of segment for building structure
	 //MeshBase::node_iterator nd = _as->get_device()->get_mesh().nodes_begin();
     //    const MeshBase::node_iterator nd_end = _as->get_device()->get_mesh().nodes_end();
	 elem = (*_structure_elements.begin());
	  nd = elem->get_node(1);

	 double edge_min, edge_max;

	 edge_min = (*nd)(0); edge_max = (*nd)(0);
	 nd++;

	 for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){
	//for ( ; nd != nd_end; nd++ ){

		 elem = *it;

		 for (unsigned int i = 0; i < elem->n_nodes(); i++){
			 nd = elem->get_node(i);
		 if ( (*nd)(0) < edge_min ) edge_min = (*nd)(0);
		if ( (*nd)(0) > edge_max ) edge_max = (*nd)(0);
		 }
        };

	_local_origin(1) = edge_min; _local_origin(2) = 0.0; _local_origin(3) = 0.0;

	std::cout << "local origin " << _local_origin << std::endl;

	double l1 = (fabs(edge_max - edge_min)) * scale;

	//Minimum periodic direction is considered along y and z axis, but eventually other lenghts can be
	//specified by user in input (conventional cells along these direction are assured also in this case!!)
	double l2 = _as->get_options().get_option("y_lenght", 0.0);
	double l3 = _as->get_options().get_option("z_lenght", 0.0);

	make_supercell( l1, l2, l3);

	//print_basis(_super_basis, "supercell.xyz");

	//std::cout << "Period is " << _period << std::endl;

 }




//void AtomisticGenerator1D::passivate(){

//  std::vector<Atom> periodic_basis;
//  unsigned int i,j,k;
//  Atom tmp;
//  Tensor1 y_period, z_period;

//  y_period(1) = _period(1,2); y_period(2) = _period(2,2); y_period(3) = _period(3,2);
//  z_period(1) = _period(1,3); z_period(2) = _period(2,3); z_period(3) = _period(3,3);
//
//  for (i = 0; i < _structure_basis.size(); i++){
//
//    tmp = _structure_basis[i];
//    tmp.set_flag(1);
//    periodic_basis.push_back(tmp);
//
//    for (j = -1; j < 2; j++){
//      for (k = -1; k < 2; k++){
//
//   if ( (j == 0)&&(k == 0) ) continue;
//
//    tmp=_structure_basis[i];
//    tmp.set_position(tmp.get_position() + j * y_period + k * z_period );
//    //tmp.position = tmp.position + j * y_period + k * z_period;
//    periodic_basis.push_back(tmp);
//
//  }
//    }
//
//  }
//
//  //print_basis(periodic_basis, "unpassivated_periodic.xyz");
//  //bond_map_periodic = bond_map_gen(periodic_basis);
//  if (_bondmapobject != NULL) delete _bondmapobject;
//  passivate_cluster(periodic_basis);
//
//  _structure_basis.clear();
//
//  for (i = 0; i < periodic_basis.size(); i++){
//    if (periodic_basis[i].get_flag() == 1) _structure_basis.push_back(periodic_basis[i]);
//
//
//}
//  //print_basis(periodic_basis, "passivated_periodic.xyz");
//  //print_basis(_structure_basis, "passivated_basis.xyz");
//
//}

// void AtomisticGenerator1D::passivate(){
//
//passivate_cluster(_structure_basis);
//}

//  Tensor1 u, u1, u2, r1, r2, r3, O;
//  double R1, R2;
//  const double sq3 = sqrt(3.0);
//  const double d = 0.64;
//  unsigned int i;
//  std::vector<Atom> hydrogens;
//  Atom tmp;
//  double sin109, cos109;
//  unsigned int ** bond_map;
//  std::vector<Atom> &basis = _structure_basis;
//
//  sin109 = sin ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );
//  cos109 = cos ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );
//
//
//
//  if (_bondmapobject == NULL) {
//    bond_map_gen(basis);
//  }
//
//  bond_map = _bondmapobject->get_bond_map();
//
//  for (i = 0; i < basis.size(); i++){
//
//    if ( bond_map[i][8] == 3 ){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//      r2 = basis[ bond_map[i][1] ].get_position() - O;
//      r3 = basis[ bond_map[i][2] ].get_position() - O;
//
//      u1 = r1 + r2 + r3;
//      R1 = norm(u1);
//      u1 = -d * (u1/R1);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      tmp.set_position ( O + u1 );
//      hydrogens.push_back(tmp);
//    }
//
//    else if  ( bond_map[i][8] == 2 ){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//      r2 = basis[ bond_map[i][1] ].get_position() - O;
//
//      u1 = r1 + r2;
//      u2 = vectorProduct(r1, r2);
//      R1 = norm(u1); R2 = norm(u2);
//
//      u = - (u1 / R1) - sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//
//      u = - (u1 / R1) + sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//    }
//
//    else if (bond_map[i][8] == 1){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//
//      Tensor2Gen vect_rot(0);
//      vect_rot(1,1) = cos109; vect_rot(1,2) = sin109; vect_rot(1,3) = 0.0;
//      vect_rot(2,1) = sin109; vect_rot(2,2) = cos109; vect_rot(2,3) = 0.0;
//      vect_rot(3,1) = 0.0; vect_rot(3,2) = 0.0; vect_rot(3,3) = 1.0;
//
//      u = vect_rot * r1;
//      R1 = norm(r1);
//      r2 = u * R1;
//      u = d * (u / norm(u));
//
//      tmp.set_specie("H");
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//      //u = (-r1)/norm(r1); u(1) = u(1) + 0.1;
//      //////////////////////////////
//
//
//
//      tmp.set_position ( O + u );
//      tmp.set_flag ( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      hydrogens.push_back(tmp);
//
//      u1 = r1 + r2;
//      u2 = vectorProduct(r1, r2);
//      R1 = norm(u1); R2 = norm(u2);
//
//      u = - (u1 / R1) - sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//           // u = (-r1)/norm(r1); u(2) = u(2) + 0.1;
//            //////////////////////////////
//
//      tmp.set_position ( O + u );
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      hydrogens.push_back(tmp);
//
//      u = - (u1 / R1) + sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//                // u = (-r1)/norm(r1); u(2) = u(2) - 0.1; u(1)=u(1) - 0.1;
//                 //////////////////////////////
//
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//    }
//
//    else if (bond_map[i][8] != 4) {std::cout << "Warning! atom " << i
//      << " is bonded to " << bond_map[i][8] << " atoms" << std::endl;}
//
//  }
//  for (i = 0; i < hydrogens.size(); i++) {basis.push_back( hydrogens[i] );}




 //}
