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

	ax = 0.0; ay = 0.0; az = 0.0;
}


AtomisticGenerator2D::~AtomisticGenerator2D(void) {};


void 
 AtomisticGenerator2D::build()
 {
   double min_x, min_y, max_x, max_y;

//Common building operations
	 make_conv_cell();
	 make_conv_basis();

	 unsigned int dimension = _as->get_device()->get_mesh().mesh_dimension(); 
		 
	 //Check edges of segment for building structure
	 MeshBase::node_iterator nd = _as->get_device()->get_mesh().nodes_begin();	
         const MeshBase::node_iterator nd_end = _as->get_device()->get_mesh().nodes_end();


	 min_x = (**nd)(0); max_x = (**nd)(0); min_y = (**nd)(1); max_y = (**nd)(1);
	 nd++;

	 //double edge_min, edge_max;
                   //edge_min = (**nd)(0); edge_max = (**nd)(0);
	 nd++;
	 
	for ( ; nd != nd_end; nd++ ){
		if ( (**nd)(0) < min_x ) min_x = (**nd)(0);                              
		if ( (**nd)(0) > max_x ) max_x = (**nd)(0);
		if ( (**nd)(1) < min_y ) min_y = (**nd)(1);                              
		if ( (**nd)(1) > max_y ) max_y = (**nd)(1);	
        };
	
	//_local_origin(1) = edge_min; _local_origin(2) = 0.0; _local_origin(3) = 0.0;
	_local_origin(1) = min_x; _local_origin(2) = min_y; _local_origin(3) = 0.0;

	double l1 = (max_x - min_x) * scale; 
	double l2 = (max_y - min_y) * scale;

	//Minimum periodic direction is considered along z axis, but eventually other lenghts can be 
	//specified by user in input (conventional cells along these direction are assured also in this case!!)
	//double l2 = _as->get_options().get_option("y_lenght", 0.0);
	double l3 = _as->get_options().get_option("z_lenght", 0.0);
   
	make_supercell( l1, l2, l3, true, true);

	//print_basis(_super_basis, "supercell.xyz");
	
	//std::cout << "Period is " << _period << std::endl;
    
 }


void AtomisticGenerator2D::passivate(void){};
