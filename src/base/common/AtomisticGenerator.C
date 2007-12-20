#include "AtomisticGenerator.h"
#include "AtomisticStructure.h"
#include "AtomisticGenerator1D.h"


const double AtomisticGenerator::tol = 1e-3;

const double AtomisticGenerator::scale = 1e2;


// AtomisticGenerator::AtomisticGenerator(AtomisticStructure* const as)
// {
//   _as = as;
// 	_rotation(1,1) = 1.0; _rotation(1,2) = 0.0; _rotation(1,3) = 0.0; _rotation(2,1) = 0.0; _rotation(2,2) = 1.0;
// 	_rotation(2,3) = 0.0; _rotation(3,1) = 0.0; _rotation(3,2) = 0.0; _rotation(3,3) = 1.0;
// }


 AtomisticGenerator* 
 AtomisticGenerator::create(AtomisticStructure* const as, unsigned int dimension)
 { 
   AtomisticGenerator* ag =  NULL;

   if (dimension == 1)  ag = AtomisticGenerator1D::create(as);
   return ag;
 }


void AtomisticGenerator::print_basis(std::vector<Atom> &basis, const std::string filename){

  std::ofstream output_file;

  std::vector<Atom>::iterator basis_iterator = basis.begin();

  output_file.open(filename.c_str());
  output_file << basis.size() << std::endl << std::endl;
  std::cout << "Began printing file " << std::endl;
  do{

    output_file << std::setw(2) << (*basis_iterator).specie 
		<< std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).position(1)) 
                                    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).position(2)) 
		<< std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).position(3)) << "\n"; 

    basis_iterator++;

  }while(basis_iterator != basis.end());
  std::cout << "Finished, I'm closing file" << std::endl;
  output_file.close();
};


void 
AtomisticGenerator::do_init()
{

  std::cout << "Beginning AtomisticGenerator::do_init()" << std::endl;
  
	
// Set material informations 
//-----------------------------------------------------------------------------------------
  std::string structure;
  structure = "none";

 std::string material;
 Material* region_material = NULL;
	if ( _as->get_options().find_option("material")) {
    material = _as->get_options().get_option("material", "Si");
  }
  else
    {
      if (!(_as->get_options().find_option("reference_region"))){
	    std::cerr << "No material could be set " << std::endl;}
     	  std::vector<ID> ids;
	  std::string ref_region;
	  ref_region = _as->get_options().get_option("reference_region", "None");
	  _as->get_device()->get_region_ids(ref_region, ids);
	  region_material = _as->get_device()->get_material(ids[0]);
	  material = region_material->get_name();
                    structure =  region_material->get_structure();
    }

  if (_as->get_options().find_option("structure")){
	_as->get_options().get_option("structure", "zb");
      }

  //Set appropriate options for corresponding materials and structures
  if ( material.compare("Si") == 0 ) 
    {std::cout << "Setting silicon options " << std::endl;
      ax = 5.4307; ay = ax; az = ax; set_lattice_type("fcc"); set_crystal_basis("zincblende", "A", "A");}
  else if (material.compare("GaAs") == 0)
    {ax = 5.6535;  ay = ax; az = ax; set_lattice_type("fcc"); set_crystal_basis("zincblende", "A", "B");}
  else if ( material.compare("Diamond") == 0 ) 
    {ax =3.56685; ay = ax; az = ax; set_lattice_type("fcc"); set_crystal_basis("diamond", "A");}
 else if ( material.compare("GaN") == 0 ) 
   {ax =3.190; ay = ax; az = 5.190; set_lattice_type("hexagonal"); set_crystal_basis("wurtzite", "A", "B", 0.0);}
    //...to be completed...


    
      
    
    
    
//Set Growth direction informations

    Tensor2Gen miller(1);
    miller(2,2) = 0.0; 
    miller(3,3) = 0.0;
    
    std::vector<int> growth_direction;

    if (! ( (region_material->get_options().find_option("x_growth_direction")) || (_as->get_options().find_option("x_growth_direction")) ) )
      std::cerr << "Warning, no x_growth_direction is set for atomistic structure " << _as->get_name() << " Setting (1,0,0) as default " << std::endl;

    if (region_material != NULL){
	if (region_material->get_options().find_option("x_growth_direction")){
region_material->get_options().get_option("x_growth_direction", growth_direction);
		miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
	}
	if (region_material->get_options().find_option("y_growth_direction")){
region_material->get_options().get_option("y_growth_direction", growth_direction);
		miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
	}
	if (region_material->get_options().find_option("z_growth_direction")){
region_material->get_options().get_option("z_growth_direction", growth_direction);
		miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
	}
    }
    
    //If Miller indexes specified in Atomistic options, take them
  	if (_as->get_options().find_option("x_growth_direction")){
region_material->get_options().get_option("x_growth_direction", growth_direction);
		miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
	}
	if (_as->get_options().find_option("y_growth_direction")){
region_material->get_options().get_option("y_growth_direction", growth_direction);
		miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
	}
	if (_as->get_options().find_option("z_growth_direction")){
region_material->get_options().get_option("z_growth_direction", growth_direction);
		miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
	}  
    
	set_prim_miller(miller);

    //-------------------------------------------------------------------------------------------
    


    // Set the vectore of elements covered by structure, useful for change specie and cut
   MeshBase::element_iterator el = _as->get_device()->get_mesh().elements_begin();	
   const MeshBase::element_iterator el_end = _as->get_device()->get_mesh().elements_end();

    for ( ; el != el_end; el++)
      {           Elem* elem = *el;
	if (_as->get_IDset().find( elem->subdomain_id() ) != _as->get_IDset().end() ) _structure_elements.push_back(elem); 
      }
   


//Build up structure with proper options
//----------------------------------------------------------------------------------------------
  build();
	

  std::string preserve;
  preserve = _as->get_options().get_option("preserve", "none");
 change_specie(preserve);


 std::string passivation;
 passivation = _as->get_options().get_option("passivation", "no");

 if (passivation.compare("yes") == 0){passivate();}


//----------------------------------------------------------------------------------------------
	
	
//Pass data to AtomisticStructure	
//--------------------------------------------------------------------------------------------------

Atom tmp_atom;
 _as->set_structure_atoms(_structure_basis);
for (int i = 0; i < 3 ; i++){
	for (int j = 0; j < 3 ; j++){
_as->_periodicity_vectors[i][j] = _period(i+1,j+1);
	}
}
_as->N_atoms = _structure_basis.size();

 std::set<std::string> atom_types;
 for (int i = 0; i < _structure_basis.size(); i++){
   atom_types.insert(_structure_basis[i].specie);
 }
 _as->N_types = atom_types.size();
 _as->_atom_types.clear();
 for (std::set<std::string>::iterator types = atom_types.begin(); types != atom_types.end(); types++){
_as->_atom_types.push_back(*types);
 }


#ifdef DEBUG	
  std::cout << "Ending AtomisticGenerator::do_init() " << std::endl;
#endif

 };



void 
AtomisticGenerator::change_specie(std::string preserve){
 
  std::set<ID> IDs = _as->get_IDset();

  std::map<std::string , std::string> assign;

  _structure_basis.clear();

  for (std::set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); reg++){

    std::string material;
    Material* material_pointer = NULL;

    material_pointer = _as->get_device()->get_material( (*reg) );
    material = material_pointer->get_name();
    assign.clear();

    //If some doping is present, a strategy must be studied and implemented here,
    // as we can choose arbitrarily species name (e.g. calling one doping Silicon Si1 and another one Si2)
    if ( material.compare("Si") == 0 ) {assign["A"] = "Si";}
    else if ( material.compare("GaN") == 0) {assign["A"] = "Ga"; assign["B"] = "N";}
    else if ( material.compare("GaAs") == 0) {assign["A"] = "Ga"; assign["B"] = "As";}
    else if ( material.compare("AlAs") == 0) {assign["A"] = "Al"; assign["B"] = "As";}
    else if ( material.compare("Diamond") == 0 ) {assign["A"] = "C";}
 
    const std::map<std::string, std::string>::iterator assign_last = assign.end();
 
    std::cout << "Assign[A] " << assign["A"] << std::endl; 
    std::cout << "Assign[B] " << assign["B"] << std::endl; 
    //Cycle upon all atoms and change specie according to assign map
    Point p(0.0, 0.0, 0.0);


    //Different strategies if preserving conventional cell or preserving basis are needed
    if (preserve.compare("none") == 0)
      {

	for ( std::vector<Atom>::iterator atom = _super_basis.begin(); atom != _super_basis.end(); atom++){
	  for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){
	    p(0) = (*atom).position(1) / scale;
	    Elem* elem = *it;
	    if ( elem->subdomain_id() == *reg){
	      if ( elem->contains_point(p) ) {
		if ( assign.find( (*atom).specie ) != assign_last ){
		  std::string tmp =  assign[(*atom).specie];
		  (*atom).specie = tmp;
		  (*atom).id = *reg;
		  _structure_basis.push_back(*atom); } 
	      }
	    }
	  }

	}
      }


    if (preserve.compare("lattice") == 0)
      {

	Atom tmp_atom;
	for ( std::vector<Tensor1>::iterator lattice = _super_lattice.begin(); lattice != _super_lattice.end(); lattice++){
	  for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){
	    p(0) = (*lattice)(1) / scale;
	    Elem* elem = *it;
	    if ( elem->subdomain_id() == *reg){
	      if ( elem->contains_point(p) ) {

		for ( std::vector<Atom>::iterator atom = _crystal_basis.begin(); atom != _crystal_basis.end(); atom++){

		  //tmp_atom.specie = (*atom).specie;
		  //		    std::cout << "specie" << tmp_atom.specie <<std::endl;
		  tmp_atom.position=_local_origin +  (*lattice) + _rotation*_prim_vec*(*atom).position;
		  tmp_atom.id = *reg;
		  if ( assign.find( (*atom).specie ) != assign_last ){
		    std::string tmp =  assign[(*atom).specie];
		    tmp_atom.specie = tmp;
		    _structure_basis.push_back(tmp_atom); } 

		}
	      }
	    }
	  }

	}
      }

    if (preserve.compare("conventional") == 0)
      {

	Atom tmp_atom;           
	for ( std::vector<Tensor1>::iterator conv = _super_conv.begin(); conv != _super_conv.end(); conv++){
	   for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++){
	    p(0) = (*conv)(1) / scale;
	    Elem* elem = *it;
	    if ( elem->subdomain_id() == *reg){
	      if ( elem->contains_point(p) ) {
		for ( std::vector<Tensor1>::iterator conv_lattice_basis_it = _conv_lattice_basis.begin(); conv_lattice_basis_it != _conv_lattice_basis.end(); conv_lattice_basis_it++){
		  for ( std::vector<Atom>::iterator atom = _crystal_basis.begin(); atom != _crystal_basis.end(); atom++){
		    tmp_atom.position=_local_origin +  (*conv) + (*conv_lattice_basis_it) + _rotation*_prim_vec*(*atom).position;
		    tmp_atom.id = *reg;
		    if ( assign.find( (*atom).specie ) != assign_last ){
		      std::string tmp =  assign[(*atom).specie];
		      tmp_atom.specie = tmp;
		      _structure_basis.push_back(tmp_atom); } 
		  }
		}
	      }
	    }

	  }
	}
      }

      
  }
}





void AtomisticGenerator::set_cutoff()
{
  _cutoff["Si"] = 1.81;
};



 
 void AtomisticGenerator::make_supercell(double l1, double l2, double l3, bool preserve_basis, bool preserve_conv){

  //Build a supercell, defined by the lenght of conventional growth cell vectors
  std::vector<Tensor1>::iterator conv_iterator;
  std::vector<Atom>::iterator basis_iterator;
  int i,j,l;
  int n1,n2,n3;
  double conv_l1, conv_l2, conv_l3;
  Atom basis_atom;
  Tensor1 lattice_point;
  Tensor2Gen supercell_vect,inv_supercell_vect;
  Tensor1 tmp_check, tmp_conv;
  bool check_boundary, check_boundary2;

 
  //Check values. l1,l2,l3 cannot be unwisely large (no more than (1um)^3)
  assert((l1*l2*l3) < 1e+12);

  //Find lenght of conventional cell sides
  conv_l1 = sqrt(_conv_vect(1,1) * _conv_vect(1,1) + _conv_vect(2,1) * _conv_vect(2,1) + _conv_vect(3,1) * _conv_vect(3,1));
  conv_l2 = sqrt(_conv_vect(1,2) * _conv_vect(1,2) + _conv_vect(2,2) * _conv_vect(2,2) + _conv_vect(3,2) * _conv_vect(3,2));
  conv_l3 = sqrt(_conv_vect(1,3) * _conv_vect(1,3) + _conv_vect(2,3) * _conv_vect(2,3) + _conv_vect(3,3) * _conv_vect(3,3));

  n1 = int(floor(l1 / conv_l1)); n2 = int(floor(l2 / conv_l2)); n3 = int(floor(l3 / conv_l3));

  if (preserve_conv) {

    l1 = (n1 + 1) * conv_l1; l2 = (n2 +1) * conv_l2; l3 = (n3 + 1) * conv_l3;

}


  //Set supercell periodical vectors
  Tensor2Gen lmat(0);

  // Periodicity along x direction is set to a big value (double of structure lenght) (non periodic along x)
  lmat(1,1) = (n1 + 1) * 2; 
  lmat(2,2) = (n2 + 1); lmat(3,3) = (n3 +1);

  _period = _conv_vect * lmat;

  //Define vectors with same direction of conventional cell vectors, but with size specifed by l1,l2,l3
  supercell_vect(1,1) = _conv_vect(1,1) * (l1 / conv_l1); supercell_vect(2,1) = _conv_vect(2,1) * (l1 / conv_l1); supercell_vect(3,1) = _conv_vect(3,1) * (l1 / conv_l1);
  supercell_vect(1,2) = _conv_vect(1,2) * (l2 / conv_l2); supercell_vect(2,2) = _conv_vect(2,2) * (l2 / conv_l2); supercell_vect(3,2) = _conv_vect(3,2) * (l2 / conv_l2);
  supercell_vect(1,3) = _conv_vect(1,3) * (l3 / conv_l3); supercell_vect(2,3) = _conv_vect(2,3) * (l3 / conv_l3); supercell_vect(3,3) = _conv_vect(3,3) * (l3 / conv_l3);
  inv_supercell_vect = inv(supercell_vect);

  for (i = -1 ; i <= n1 + 1 ; i++){
    for (j = -1 ; j <= n2 + 1 ; j++){
      for (l = -1 ; l <= n3 + 1 ; l++){

	conv_iterator = _conv_lattice_basis.begin();

	//Fill conventional edges basis (super_conv)
	if ( (i != -1)&&(i <= n1 )&&(j != -1)&&(j <= n2)&& (l != -1)&&(l <= n3) ){
	tmp_conv(1) = (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
	tmp_conv(2) = (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
	tmp_conv(3) = (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
	_super_conv.push_back(tmp_conv);}

	do{ 
	  //Assign lattice point position 
	  lattice_point(1) = (*conv_iterator)(1) + (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
	  lattice_point(2) = (*conv_iterator)(2) + (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
	  lattice_point(3) = (*conv_iterator)(3) + (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));

	  if (preserve_basis){

	    //Check if lattice point is inside bonduary
	     if ((i >= n1) || (i <= 0) || (j >= n2) || (j<= 0) || (l >= n3) || (l <= 0)) {
	      tmp_check = inv_supercell_vect * lattice_point;
	      check_boundary = ((tmp_check(1) >= -tol) && (tmp_check(1) <(1.0 - tol))) &&
		((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 - tol))) &&
		((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 - tol)));
	     } 

	    else (check_boundary = 1);
	   
	    if (check_boundary){
	      //Put lattice point into supercell lattice points array
	      _super_lattice.push_back(lattice_point);

	      basis_iterator=_crystal_basis.begin();

	      do{
		basis_atom = (*basis_iterator);
                                    basis_atom.position = lattice_point+_rotation*_prim_vec*(*basis_iterator).position;
		
		_super_basis.push_back(basis_atom);
		basis_iterator++;

	      }while(basis_iterator != _crystal_basis.end());
	    }
	  }

	  else{

                      
	    //Put lattice point into supercell lattice points array
	    _super_lattice.push_back(lattice_point);

	    basis_iterator=_crystal_basis.begin();

	    do{
	      basis_atom = (*basis_iterator);
	      basis_atom.position = lattice_point + _rotation*_prim_vec*(*basis_iterator).position;
		
	      //Check if basis atom is inside bonduary when preserve_basis is off
	      if ((i >= n1) || (i <= 0) || (j >= n2) || (j<= 0) || (l >= n3) || (l <= 0)) {
		tmp_check = inv_supercell_vect * basis_atom.position;
		check_boundary2 = ((tmp_check(1) >= -tol) && (tmp_check(1) <(1.0 + tol))) &&
		  ((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 + tol))) &&
		  ((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 + tol)));
		if (check_boundary2) _super_basis.push_back(basis_atom);
	      } 
	      else _super_basis.push_back(basis_atom);
	      basis_iterator++;
	    }while(basis_iterator != _crystal_basis.end());
	  }
              

	  conv_iterator++;
	}while(conv_iterator != _conv_lattice_basis.end());
     
      };
    };
  };
};


void 
AtomisticGenerator::set_lattice_type(const std::string lattice_name)
{

  Tensor2Gen prim_vec_dir;

  //Set the lattice type. It defines the primitive vectors
  _lattice_type=lattice_name;
    
  if (_lattice_type.compare("cubic") == 0) {

    assert((ax == ay) && (ay == az));

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0; prim_vec_dir(3,1) = 0;
    prim_vec_dir(1,2) = 0; prim_vec_dir(2,2) = 1; prim_vec_dir(3,2) = 0;
    prim_vec_dir(1,3) = 0; prim_vec_dir(2,3) = 0; prim_vec_dir(3,3) = 1;

    _prim_vec = prim_vec_dir * ax;

  } 

  else if (_lattice_type.compare("bcc") == 0) {

    assert((ax == ay) && (ay == az));

    prim_vec_dir(1,1) = -0.5; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = -0.5; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = -0.5;

    _prim_vec = prim_vec_dir * ax;

  } 

  else if (_lattice_type.compare("fcc") == 0) {

    assert((ax == ay) && (ay == az));

    prim_vec_dir(1,1) = 0.0; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = 0.0; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.0;

    _prim_vec = prim_vec_dir * ax;

  } 

  else if (_lattice_type.compare("hexagonal") == 0) {

    assert(ax == ay);

    prim_vec_dir(1,1) = 0.5; prim_vec_dir(2,1) = -sqrt(3.0) / 2.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = sqrt(3.0) / 2.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.0; prim_vec_dir(2,3) = 0.0; prim_vec_dir(3,3) = 1.0;

    _prim_vec(1,1) = prim_vec_dir(1,1) * ax; _prim_vec(2,1) = prim_vec_dir(2,1) * ax; 
    _prim_vec(1,2) = prim_vec_dir(1,2) * ax; _prim_vec(2,2) = prim_vec_dir(2,2) * ax; 
    _prim_vec(1,3) = 0.0; _prim_vec(2,3) = 0.0; _prim_vec(3,3) = prim_vec_dir(3,3) * az;

  } 

  else if (_lattice_type.compare("anatase") == 0) {

    assert(ax == ay);

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.5;

    _prim_vec(1,1) = prim_vec_dir(1,1) * ax; 
    _prim_vec(2,2) = prim_vec_dir(2,2) * ax; 
    _prim_vec(1,3) = prim_vec_dir(1,3) * ax; _prim_vec(2,3) = prim_vec_dir(2,3) * ax; _prim_vec(3,3) = prim_vec_dir(3,3) * az;

}

  else std::cout << "Lattice type " << _lattice_type << " doesn't exist" << std::endl;
  
};


void AtomisticGenerator::set_crystal_basis(const std::string basis_name, const std::string specie1, const std::string specie2, double u)
{
  //Set the atomic basis.
  //Vectors are defined in primitive vectors basis.
  _basis_type = basis_name;
  Atom tmp;

  assert(_crystal_basis.empty());

  if (_basis_type.compare("cubic") == 0){
    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
     tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("diamond") == 0){
 
    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
     tmp.position(1) = -1.0 / 8.0; tmp.position(2) = -1.0/8.0; tmp.position(3) = -1.0/8.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1)=1.0/8.0; tmp.position(2)=1.0/8.0; tmp.position(3)=1.0/8.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("zincblende") == 0){

    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);
 
    if (specie1.compare("not_specified") == 0) tmp.specie = "B";
    else tmp.specie = specie2;
    tmp.position(1) = 1.0 / 4.0; tmp.position(2) = 1.0 / 4.0; tmp.position(3) = 1.0 / 4.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("cristobalite") == 0){

    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
    tmp.position(1) = 1.0 / 8.0; tmp.position(2) = 1.0 / 8.0; tmp.position(3) = 1.0 / 8.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = -1.0 / 8.0; tmp.position(2) = -1.0 / 8.0; tmp.position(3) = -1.0 / 8.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "B";
    else tmp.specie = specie2;
    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 1.0 / 2.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 0.0; tmp.position(2) = 1.0 / 2.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 1.0 / 2.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("NaCl") == 0){

    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "B";
    else tmp.specie = specie2;
    tmp.position(1) = 1.0 / 2.0; tmp.position(2) =1.0 / 2.0; tmp.position(3) = 1.0 / 2.0;
    _crystal_basis.push_back(tmp);

  }

  else if(_basis_type.compare("wurtzite") == 0){

    if (u == 0.0) u = 3.0 / 8.0;
    if (specie1.compare("not_specified") == 0) tmp.specie = "A";
    else tmp.specie = specie1;
    tmp.position(1) = 1.0 / 3.0; tmp.position(2) = 2.0 / 3.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 2.0 / 3.0; tmp.position(2) = 1.0 / 3.0; tmp.position(3) = 1.0 / 2.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "B";
    else tmp.specie = specie2;
     tmp.position(1) = 1.0 / 3.0; tmp.position(2) = 2.0 / 3.0; tmp.position(3) = u;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 2.0 / 3.0; tmp.position(2) = 1.0 / 3.0; tmp.position(3) = 1.0 / 2.0 + u;
    _crystal_basis.push_back(tmp); 
  }

  else if (_basis_type.compare("TiO2") == 0) {

   if (u == 0.0) u = 0.25;
   if (specie1.compare("not_specified") == 0) tmp.specie = "A";
   else tmp.specie = specie1;
   tmp.position(1) = -1.0 / 8.0; tmp.position(2) = 5.0 / 8.0; tmp.position(3) = 1.0 / 4.0;
   _crystal_basis.push_back(tmp);

   tmp.position(1) = 1.0 / 8.0; tmp.position(2) = 3.0 / 8.0; tmp.position (3) = 3.0 / 4.0;
   _crystal_basis.push_back(tmp);

   if (specie1.compare("not_specified") == 0) tmp.specie = "B";
   else tmp.specie = specie1;
   tmp.position(1) = -u; tmp.position(2) = 0.25 - u; tmp.position(3) = 2.0 * u;
   _crystal_basis.push_back(tmp);

  tmp.position(1) = - u  -  0.25; tmp.position(2) = 0.5 - u; tmp.position(3) = 0.5 + 2.0 * u;
   _crystal_basis.push_back(tmp);

 tmp.position(1) = + u; tmp.position(2) = - 0.25 - u; tmp.position(3) = 1.0 - 2.0 * u;
   _crystal_basis.push_back(tmp);

 tmp.position(1) = + 0.25 + u; tmp.position(2) =  0.5 + u; tmp.position(3) = 0.5 - 2.0 * u;
   _crystal_basis.push_back(tmp);

}

  else std::cout << "Crystal atom basis specified does not exist" << std::endl;

assert(~(_crystal_basis.empty()));

#ifdef DEBUG
 std::cerr << "done" << std::endl;
#endif
};





void AtomisticGenerator::set_prim_miller(Tensor2Gen cut_planes)
{
	//Reduce cut planes if they're not in minimal integer form
 
	cut_planes = reduce_vector(cut_planes);
 
  //Transform miller indexes in primitive reciprocal vectors basis 
  //Needed when conventional cell differs from primitive cell
  Tensor2Gen prim_miller_basis = reciprocal(_prim_vec);

  if (_lattice_type.compare("cubic") == 0){
    //In cubic lattice conventional cell equals to unit cell
    _prim_miller = cut_planes / (ax);
  }

  else if (_lattice_type.compare("bcc") == 0){
    _prim_miller = (inv(prim_miller_basis) * cut_planes / (ax));
    scale_to_int(_prim_miller);
  }

  else if (_lattice_type.compare("fcc") == 0){
    _prim_miller=inv(prim_miller_basis) * cut_planes / (ax);
    scale_to_int(_prim_miller);
  }

  else if (_lattice_type.compare("hexagonal") == 0){
    _prim_miller(1,1) = cut_planes(1,1) / ax; _prim_miller(2,1) = cut_planes(2,1) / ay; _prim_miller(3,1) = cut_planes(3,1) / az;
    _prim_miller(1,2) = cut_planes(1,2) / ax; _prim_miller(2,2) = cut_planes(2,2) / ay; _prim_miller(3,2) = cut_planes(3,2) / az;
    _prim_miller(1,3) = cut_planes(1,3) / ax; _prim_miller(2,3) = cut_planes(2,3) / ay; _prim_miller(3,3)=cut_planes(3,3) / az;
    scale_to_int(_prim_miller);
  }

  else {
    _prim_miller(1,1) = cut_planes(1,1) / ax; _prim_miller(2,1) = cut_planes(2,1) / ay; _prim_miller(3,1) = cut_planes(3,1) / az;
    _prim_miller(1,2) = cut_planes(1,2) / ax; _prim_miller(2,2) = cut_planes(2,2) / ay; _prim_miller(3,2) = cut_planes(3,2) / az;
    _prim_miller(1,3) = cut_planes(1,3) / ax; _prim_miller(2,3) = cut_planes(2,3) / ay; _prim_miller(3,3) = cut_planes(3,3) / az;
    scale_to_int(_prim_miller);
};

 
};








void AtomisticGenerator::make_conv_cell()
{
  //Calculate conventional cell vectors in the directions given by cut planes (conventional growth cell)
  Tensor1 m1,m2,m3,select_vect(0);
  Tensor1 conv1, conv2, conv3;
  int i;
  Tensor2Gen prim_miller_basis;

  prim_miller_basis = reciprocal(_prim_vec);
  _conv_prim = inv(_prim_vec) * prim_miller_basis * _prim_miller;
  scale_to_int(_conv_prim);
  _conv_vect = _prim_vec * _conv_prim;

  for (i = 1; i <= 3; i++) {conv1(i) = _conv_vect(i, 1); conv2(i) = _conv_vect(i, 2); conv3(i) = _conv_vect(i, 3);};

  //At least x growth index must be specified (growth direction in 1D structure)
  assert (norm(conv1) > tol);

  // If no other indexes are specified, set y growth direction by default (orthogonal to x index)
  if ((norm(conv2) < tol)&&(norm(conv3) < tol)) {

   std::cerr << "Warning: only x growht direction is defined. Building other direction orthogonal " << std::endl;

  conv2(1) = - conv1(3); conv2(2) = 0.0; conv2(3) = conv1(1);
  conv3 =  vectorProduct(conv1, conv2);
  conv3 = conv3 / norm(conv3);
}

  //If only one index is not specified, build it by default (orthogonal)
  else if (norm(conv2) < tol){

std::cerr << "Warning: only x and z growht direction is defined. Building y direction orthogonal " << std::endl;

    conv2 =  vectorProduct(conv1, conv3);
    conv2 = conv2 / norm(conv2);
  }

  else if (norm(conv3) < tol){

std::cerr << "Warning: only x and y growht direction is defined. Building z direction orthogonal " << std::endl;

    conv3 =  vectorProduct(conv1, conv2);
    conv3 = conv3 / norm(conv3);
  };


//   //If not all miller indexes are specified, build others as a default 
//   if (norm(conv1) < tol) {conv1(1) = - conv3(3); conv1(2) = 0.0; conv1(3) = conv3(1);}
  
//   //If only z and x miller indexes are specifed, get y miller index
//   if (norm(conv2) < tol) {
//     conv2 =  vectorProduct(conv1, conv3);
//     conv2 = conv2 / norm(conv2);
//   };


  for (i = 1; i <= 3; i++) {_conv_vect(i, 1) = conv1(i); _conv_vect(i, 2) = conv2(i); _conv_vect(i, 3) = conv3(i);};
  

  _conv_prim = inv(_prim_vec) * _conv_vect;
  scale_to_int(_conv_prim);
  _conv_vect = _prim_vec * _conv_prim;

  // Calculate distance between equivalent planes for every cut plane
  //prim_miller = inv(prim_miller_basis) * prim_vec * conv_prim;
  //for (i = 1; i <= 3; i++) {m1(i) = prim_miller(i, 1); m2(i) = prim_miller(i, 2); m3(i) = prim_miller(i, 3);};
  //planar_distance[0] = 1 / norm(prim_miller_basis * m1);
  //planar_distance[1] =1 / norm(prim_miller_basis * m2);
  //planar_distance[2] = 1 / norm(prim_miller_basis * m3);

};


void AtomisticGenerator::make_conv_basis()
{
  //Fill the conventional growth cell with atomic basis
  int lower_1, lower_2, lower_3, upper_1, upper_2, upper_3, i;
  Tensor1 prim_position, tmp_check;
  Tensor1 tmp_position;
  Tensor1 vec_x,vec_y,vec_z;
  std::vector<Atom>::iterator basis_iterator;
  Tensor2Gen rot_tmp;

  //Make a preliminar rotation only if conventional cell vectors are orthogonal
 
  vec_x(1) = _conv_vect(1,1); vec_x(2) = _conv_vect(2,1); vec_x(3) = _conv_vect(3,1);
  vec_y(1) = _conv_vect(1,2); vec_y(2) = _conv_vect(2,2); vec_y(3) = _conv_vect(3,2);
  vec_z(1) = _conv_vect(1,3); vec_z(2) = _conv_vect(2,3); vec_z(3) = _conv_vect(3,3);


 //If z and y vectors are not orthogonal, build rotation vector y as a vector orthogonal to x
  if ( (vec_z * vec_y ) > 1e-10) {

    if ( (vec_y * vec_x) < 1e-10 ) {vec_z = vectorProduct(vec_x, vec_y);}
    else if  ( (vec_z * vec_x) < 1e-10)  {vec_y = vectorProduct(vec_x, vec_z);}
    else {std::cout << "Warning: at least x or y growth direction should be orthogonal to z growth direction" << std::endl;}

  }


  //vec_y = vectorProduct(vec_x, vec_z);
  if (((vec_x * vec_y) < 1e-10) && ((vec_x * vec_z) < 1e-10) && ((vec_y * vec_z) < 1e-10)) { 

    assert(norm(vec_x) > 1e-10);
    assert(norm(vec_y) > 1e-10);
    assert(norm(vec_z) > 1e-10);

    for ( int i = 1; i <=3; i++ )  rot_tmp(1,i) = vec_x(i)/norm(vec_x);
    for ( int i = 1; i <=3; i++ )  rot_tmp(2,i) = vec_y(i)/norm(vec_y);
    for ( int i = 1; i <=3; i++ )  rot_tmp(3,i) = vec_z(i)/norm(vec_z);
  }
  else {
    rot_tmp(1,1) = 1.0; rot_tmp(1,2) = 0.0; rot_tmp(1,3) = 0.0;
    rot_tmp(2,1) = 0.0; rot_tmp(2,2) = 1.0; rot_tmp(2,3) = 0.0;
    rot_tmp(3,1) = 0.0; rot_tmp(3,2) = 0.0; rot_tmp(3,3) = 1.0;
   std::cout << "Warning: no rotation has been done, orientation of atoms in space depends on primitive vectors definitions" << std::endl;
  }
 
  //Also a user-defined rotation is allowed, which put the axes in a direction different
  //from canonical basis
  for (i = 1; i <= 3; i++){vec_x(i) = _rotation(i,1); vec_y(i) = _rotation(i,2); vec_z(i) = _rotation(i,3);}
  vec_x =vec_x / norm(vec_x); vec_y = vec_y / norm(vec_y); vec_z = vec_z / norm(vec_z);
  for (i = 1; i <= 3; i++){_rotation(i,1) = vec_x(i); _rotation(i,2) = vec_y(i); _rotation(i,3) = vec_z(i);}
  _rotation = _rotation * rot_tmp;


  //Define a box including conventional cell
  lower_1 = int(std::min(0.0,std::min(_conv_prim(1,1),std::min(_conv_prim(1,2),_conv_prim(1,3)))));
  upper_1 = int(std::max(0.0,std::max(_conv_prim(1,1),std::max(_conv_prim(1,2),_conv_prim(1,3)))));
  lower_2 = int(std::min(0.0,std::min(_conv_prim(2,1),std::min(_conv_prim(2,2),_conv_prim(2,3)))));
  upper_2 = int(std::max(0.0,std::max(_conv_prim(2,1),std::max(_conv_prim(2,2),_conv_prim(2,3)))));
  lower_3 = int(std::min(0.0,std::min(_conv_prim(3,1),std::min(_conv_prim(3,2),_conv_prim(3,3)))));
  upper_3 = int(std::max(0.0,std::max(_conv_prim(3,1),std::max(_conv_prim(3,2),_conv_prim(3,3)))));

  for (int i = lower_1 - 1; i <= upper_1 + 1; i++){
    for (int j = lower_2 - 1; j <= upper_2 + 1; j++){
      for (int l = lower_3 - 1; l <= upper_3 + 1; l++){

	prim_position(1) = double(i); prim_position(2) = double(j); prim_position(3) = double(l);
	
	tmp_position = prim_position;

	bool check_boundary;

	  tmp_check = inv(_conv_prim) * tmp_position;

	check_boundary= ((tmp_check(1) >= -tol) && (tmp_check(1) < (1.0 - tol)))&&
	  ((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 - tol))) &&
	  ((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 - tol)));

	if (check_boundary){
	  tmp_position = _rotation * _prim_vec * prim_position;
	  _conv_lattice_basis.push_back(tmp_position);
}	 

      };
    };
  };
  _conv_vect = _rotation * _conv_vect;
};



//Bond map generation (cluster)
unsigned int**  AtomisticGenerator::bond_map_gen(std::vector<Atom> &basis){

  int i, j, n, n_atom = basis.size();
  unsigned int**  bond_map;
  double cutofftmp;

#ifdef DEBUG
  std::cout << "Building bond map" << std::endl;
#endif

  set_cutoff();
 
  //Dynamic Allocation of bond_map and 0 initialization
  bond_map = new unsigned int* [n_atom];
  for (i = 0; i < n_atom; i++)  {
    bond_map[i] = new unsigned int [9];
    for (j = 0; j < 9; j++) bond_map[i][j] = 0;
  }  

  for (i = 0; i < n_atom; i ++){
  
    for (j = i + 1; j < n_atom; j ++){

    cutofftmp = _cutoff[basis[i].specie.c_str()] + _cutoff[basis[j].specie.c_str()];

    if ( cutofftmp == 0.0 ) std::cout << "WARNING, A CUTOFF DISTANCE IS NOT DEFINED " << std::endl;    
   
      {
 
      if ( norm( basis[i].position - basis[j].position) < cutofftmp ){
	
	for (n = 0; n < 8; n++){
	  if (bond_map[i][n] == 0) {bond_map[i][n] = j; bond_map[i][8]++; break;}
	}
	
	for (n = 0; n < 8; n++){
	  if (bond_map[j][n] == 0) {bond_map[j][n] = i; bond_map[j][8]++; break;} 
	}
      }

      }

    }
  }

  return bond_map;  
};



//Hydrogen passivation. Atom basis and associated bond map are needed
void AtomisticGenerator::passivate_cluster(std::vector<Atom> &basis, unsigned int** bond_map){

  Tensor1 u, u1, u2, r1, r2, r3, O;
  double R1, R2;
  const double sq3 = sqrt(3.0);
  const double d = 1.0;
  int i;
  std::vector<Atom> hydrogens;
  Atom tmp;
  double sin109, cos109;

  sin109 = sin ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );
  cos109 = cos ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );


 if (bond_map == NULL) {
bond_map = bond_map_gen(basis);
 }

  for (i = 0; i < basis.size(); i++){
   
    if ( bond_map[i][8] == 3 ){

      O = basis[i].position;
      r1 = basis[ bond_map[i][0] ].position - O;
      r2 = basis[ bond_map[i][1] ].position - O;
      r3 = basis[ bond_map[i][2] ].position - O;

      u1 = r1 + r2 + r3;
      R1 = norm(u1);
      u1 = -d * (u1/R1);

      tmp.specie = "H";
      tmp.flag = basis[i].flag;
      tmp.position = O + u1;
      hydrogens.push_back(tmp);
    }

    else if  ( bond_map[i][8] == 2 ){

      O = basis[i].position;
      r1 = basis[ bond_map[i][0] ].position - O;
      r2 = basis[ bond_map[i][1] ].position - O;

      u1 = r1 + r2;
      u2 = vectorProduct(r1, r2);
      R1 = norm(u1); R2 = norm(u2);

      u = - (u1 / R1) - sq3 * (u2 / R2);
      u = d * (u / 2.0);

      tmp.specie = "H";
      tmp.flag = basis[i].flag;
      tmp.position = O + u;
      hydrogens.push_back(tmp);

      u = - (u1 / R1) + sq3 * (u2 / R2);
      u = d * (u / 2.0);

      tmp.specie = "H";
      tmp.flag = basis[i].flag;
      tmp.position = O + u;
      hydrogens.push_back(tmp);
    }

    else if (bond_map[i][8] == 1){

      O = basis[i].position;
      r1 = basis[ bond_map[i][0] ].position - O;

      Tensor2Gen vect_rot(0);
      vect_rot(1,1) = cos109; vect_rot(1,2) = sin109; vect_rot(1,3) = 0.0;
      vect_rot(2,1) = sin109; vect_rot(2,2) = cos109; vect_rot(2,3) = 0.0;
      vect_rot(3,1) = 0.0; vect_rot(3,2) = 0.0; vect_rot(3,3) = 1.0;

      u = vect_rot * r1;
      R1 = norm(r1);
      r2 = u * R1;
      u = d * (u / norm(u));

      tmp.specie = "H";
      tmp.position = O + u;
      tmp.flag = basis[i].flag;
      hydrogens.push_back(tmp);

      u1 = r1 + r2;
      u2 = vectorProduct(r1, r2);
      R1 = norm(u1); R2 = norm(u2);

      u = - (u1 / R1) - sq3 * (u2 / R2);
      u = d * (u / 2.0);

      tmp.specie = "H";
      tmp.position = O + u;
      tmp.flag = basis[i].flag;
      hydrogens.push_back(tmp);

      u = - (u1 / R1) + sq3 * (u2 / R2);
      u = d * (u / 2.0);

      tmp.specie = "H";
      tmp.flag = basis[i].flag;
      tmp.position = O + u;
      hydrogens.push_back(tmp);
    }

    else if (bond_map[i][8] != 4) {std::cout << "Warning! An atom is bonded to 0 or more than 4 atoms" << std::endl;}

  }
  for (i = 0; i < hydrogens.size(); i++) {basis.push_back( hydrogens[i] );}



};




//Some data manipulation function useful only in this class


Tensor2Gen 
AtomisticGenerator::reciprocal(Tensor2Gen real_basis)
{
  //Build the reciprocal basis related to input 2-rank tensor
  Tensor1 a1,a2,a3;
  Tensor1 b1,b2,b3;
  Tensor1 select_vect(0);
  Tensor2Gen reciprocal;

  //Select vector a1
  select_vect(1) = 1.0; a1 = real_basis * select_vect;
  
  //Select vector a2
  select_vect(1) = 0.0; select_vect(2) = 1.0; a2 = real_basis * select_vect;
  
  //Select vector a3
  select_vect(2) = 0.0; select_vect(3) = 1.0; a3 = real_basis * select_vect;
  
  const double volume = a1 * vectorProduct(a2,a3);
  assert (volume != 0); 
  b1 = vectorProduct(a2,a3) / volume; 
  b2 = vectorProduct(a3,a1) / volume;
  b3 = vectorProduct(a1,a2) / volume;
  
  reciprocal(1,1)=b1(1); reciprocal(2,1)=b1(2);reciprocal(3,1)=b1(3);
  reciprocal(1,2)=b2(1);reciprocal(2,2)=b2(2);reciprocal(3,2)=b2(3);
  reciprocal(1,3)=b3(1);reciprocal(2,3)=b3(2);reciprocal(3,3)=b3(3);

  assert(det(reciprocal) != 0);

  return reciprocal;
};


int 
AtomisticGenerator::compare_tol(double a, double b)
{
  //Comparison routine with a tolerance defined as internal constant.
  //If absolute value of difference between a and b is minor than tolerance, 
  //a and b are considered equal
  if (std::fabs(a-b) < tol) return 1;
  else return 0;
};


int 
AtomisticGenerator::double_to_int_cast_checked(double a)
{
  //Convert a double to the nearest integer, within a certain tolerance
  int n;
  if (std::abs(std::floor(a)-a) < std::abs(std::ceil(a) - a)) n = int(std::floor(a));
  else n = int(std::ceil(a));
  assert (std::abs(double(n) - a) < tol);
  return n;
};


double 
AtomisticGenerator::double_to_int_value_checked(double a)
{
  //Gives the double number equal to the integer nearest to a, within a certain tolerance
  double b;
  b = double(double_to_int_cast_checked(a));
  return b; 
};


void 
AtomisticGenerator::double_to_int_value_checked(Tensor1& a)
{
  double tmp;;

  tmp = double_to_int_value_checked(a(1)); a(1) = tmp;
  tmp = double_to_int_value_checked(a(2)); a(2) = tmp;
  tmp = double_to_int_value_checked(a(3)); a(3) = tmp;
};


int 
AtomisticGenerator::gcd(int a, int b)
{
  //Calculate greater common denominator between integers
  //(return 0 if gcd(a,0) or gcd(0,a)

  //Added: if illegal operation gcd(a,0) performed return abs(a) (useful for reduce_vector routine)
  if ((a == 0) || (b == 0)) return std::max( std::abs(a), std::abs(b) );

  //by Derek Chandler, MEng, MIEE
  int reminder;
  do{
    reminder = a % b;
    if (reminder != 0)
      {
	a = b;
	b = reminder;
      }
  } while (reminder);
  return b;
};


Tensor1 
AtomisticGenerator::reduce_vector(Tensor1 v)
{
  //Reduce a vector of double containing integer values to its minimal form
  Tensor1 v_tmp;
  int gcd_tmp,gcd_value;

  if (norm(v) < tol) return v;

  //Find the maximimum common denominator
  gcd_tmp = gcd(double_to_int_cast_checked(v(1)),double_to_int_cast_checked(v(2)));
  gcd_value = gcd(gcd_tmp,double_to_int_cast_checked(v(3)));
  v_tmp = v / double(gcd_value);
  double_to_int_value_checked(v_tmp);
  return v_tmp;

};


void AtomisticGenerator::scale_to_int(Tensor1& a)
{
//Expand a double vector to a vector having same direction but integer values
//If vector is not in a reduced form (having integer values with gcd > 1) it's reduced
//If input is a zero vector, a zero vector is returned

//Check if a is already integer
  if ( (fabs(a(1) - round(a(1))) < tol) && (fabs(a(2) - round(a(2))) < tol) && (fabs(a(3) - round(a(3))) < tol) )
    {
      a(1) = double_to_int_value_checked(a(1));
      a(2) = double_to_int_value_checked(a(2));
      a(3) = double_to_int_value_checked(a(3));
    }

  else

    {
    
  int i = 0;
  Tensor1 a_tmp;

  do{
    i = i + 1;
    a_tmp = a * i;
   }while (  ( fabs (a_tmp(1) - round(a_tmp(1))) >= tol)  ||  (fabs(a_tmp(2) - round(a_tmp(2))) >= tol)  || (fabs(a_tmp(3) - round(a_tmp(3))) >= tol) );

  a(1) = a_tmp(1); a(2) = a_tmp(2); a(3) = a_tmp(3);

    }

  a = reduce_vector(a);

}



void 
AtomisticGenerator::scale_to_int(Tensor2Gen& a)
{
  //Same of scale_to_int with Tensor1 argument, considering the columns of a 2-rank tensor
  Tensor1 tmp;
  Tensor1 select_vect(0);

  select_vect(1) = 1.0;  tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,1) = tmp(1); a(2,1) = tmp(2); a(3,1) = tmp(3);

  select_vect(1) = 0.0; select_vect(2) = 1.0; tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,2) = tmp(1); a(2,2) = tmp(2); a(3,2) = tmp(3);

  select_vect(2) = 0.0; select_vect(3) = 1.0; tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,3) = tmp(1); a(2,3) = tmp(2); a(3,3) = tmp(3);
};


Tensor2Gen 
AtomisticGenerator::reduce_vector(Tensor2Gen a)
{
  //Same of reduce_vector with Tensor1 argument, considering the columns of a 2-rank tensor
  Tensor1 tmp1,tmp2;
  Tensor1 select_vect(0);
  Tensor2Gen b(0);

  //Select vector a1
  select_vect(1) = 1.0;  tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,1) = tmp2(1); b(2,1) = tmp2(2); b(3,1) = tmp2(3);
  
  //Select vector a2
  select_vect(1) = 0.0; select_vect(2) = 1.0; tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,2) = tmp2(1); b(2,2) = tmp2(2); b(3,2) = tmp2(3);

  //Select vector a3
  select_vect(2) = 0.0; select_vect(3) = 1.0; tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,3) = tmp2(1); b(2,3) = tmp2(2); b(3,3) = tmp2(3);
  return b;
	
};

