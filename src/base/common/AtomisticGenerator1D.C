#include "AtomisticGenerator1D.h"
#include "AtomisticStructure.h"
#include "mesh.h"



AtomisticGenerator1D* 
AtomisticGenerator1D::create(AtomisticStructure* const as)
{ 
  AtomisticGenerator1D* ag =  NULL;
  ag = new AtomisticGenerator1D(as);
  return ag;
}


AtomisticGenerator1D::AtomisticGenerator1D(AtomisticStructure* const as)
{
  _as = as;
	_rotation(1,1) = 1.0; _rotation(1,2) = 0.0; _rotation(1,3) = 0.0; _rotation(2,1) = 0.0; _rotation(2,2) = 1.0;
	_rotation(2,3) = 0.0; _rotation(3,1) = 0.0; _rotation(3,2) = 0.0; _rotation(3,3) = 1.0;

	ax = 0.0; ay = 0.0; az = 0.0;
}


 void 
 AtomisticGenerator1D::build()
 {
//Common building operations
	 make_conv_cell();
	 make_conv_basis();

	 unsigned int dimension = _as->get_device()->get_mesh().mesh_dimension(); 
		 
	 //Check edges of segment for building structure
	 MeshBase::node_iterator nd = _as->get_device()->get_mesh().nodes_begin();	
         const MeshBase::node_iterator nd_end = _as->get_device()->get_mesh().nodes_end();

	 double edge_min, edge_max;
         edge_min = (**nd)(0); edge_max = (**nd)(0);
	 nd++;
	 
	for ( ; nd != nd_end; nd++ ){
		if ( (**nd)(0) < edge_min ) edge_min = (**nd)(0);
		if ( (**nd)(0) > edge_max ) edge_max = (**nd)(0);	
        };
	
	_local_origin(1) = edge_min; _local_origin(2) = 0.0; _local_origin(3) = 0.0;
	double l1 = (edge_max - edge_min) * scale; 

	//Minimum periodic direction is considered along y and z axis, but eventually other lenghts can be 
	//specified by user in input (conventional cells along these direction are assured also in this case!!)
	double l2 = _as->get_options().get_option("y_lenght", 0.0);
	double l3 = _as->get_options().get_option("z_lenght", 0.0);
   
	make_supercell( l1, l2, l3, true, true);

	//print_basis(_super_basis, "supercell.xyz");
	
	//std::cout << "Period is " << _period << std::endl;
    
 }




void AtomisticGenerator1D::passivate(){

  std::vector<Atom> periodic_basis;
  int i,j,k;
  Atom tmp;
  unsigned int ** bond_map_periodic;
  Tensor1 y_period, z_period;

  y_period(1) = _period(1,2); y_period(2) = _period(2,2); y_period(3) = _period(3,2);
  z_period(1) = _period(1,3); z_period(2) = _period(2,3); z_period(3) = _period(3,3);
  
  for (i = 0; i < _structure_basis.size(); i++){

    tmp = _structure_basis[i];
    tmp.flag = 1;
    periodic_basis.push_back(tmp);

    for (j = -1; j < 2; j++){
      for (k = -1; k < 2; k++){

   if ( (j == 0)&&(k == 0) ) continue;

    tmp=_structure_basis[i];
    tmp.position = tmp.position + j * y_period + k * z_period;
    periodic_basis.push_back(tmp);

  }
    }

  }

  //print_basis(periodic_basis, "unpassivated_periodic.xyz");
  bond_map_periodic = bond_map_gen(periodic_basis);
  passivate_cluster(periodic_basis, bond_map_periodic);

  _structure_basis.clear();

  for (i = 0; i < periodic_basis.size(); i++){
    if (periodic_basis[i].flag == 1) _structure_basis.push_back(periodic_basis[i]);


}
  //print_basis(periodic_basis, "passivated_periodic.xyz");
  //print_basis(_structure_basis, "passivated_basis.xyz");

}
