#include "AtomisticGenerator.h"
using namespace std;

AtomAG::AtomAG(){
  flag = 0;
}

const double AtomisticGenerator::tol = 1e-14;


AtomisticGenerator::AtomisticGenerator(AtomisticStructure* const as)
{
  _as = as;
}


AtomisticGenerator* 
AtomisticGenerator::create(AtomisticStructure* const as)
{ 
  AtomisticGenerator* ag =  NULL;
  ag = new AtomisticGenerator(as);
  return ag;
}



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
	    std::cerr << "No material could be set, Silicon set by default " << std::endl;}
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
    {set_lattice_type("fcc"); set_crystal_basis("zincblende", "Si", "Si"); ax = 5.43102064; ay = ax; az = ax;}
  else if ( material.compare("Diamond") == 0 ) 
    {set_lattice_type("fcc"); set_crystal_basis("diamond", "C"); ax =3.56685; ay = ax; az = ax;}

  //...to be completed...


//Set Growth direction informations

    Tensor2Gen miller(1);
    
    std::vector<int> growth_direction;
    if (region_material != NULL){
	if (region_material->get_options().find_option("x-growth-direction")){
region_material->get_options().get_option("x-growth-direction", growth_direction);
		miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
	}
	if (region_material->get_options().find_option("y-growth-direction")){
region_material->get_options().get_option("y-growth-direction", growth_direction);
		miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
	}
	if (region_material->get_options().find_option("z-growth-direction")){
region_material->get_options().get_option("z-growth-direction", growth_direction);
		miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
	}
    }
    
    //If Miller indexes specified in Atomistic options, take them
  	if (_as->get_options().find_option("x-growth-direction")){
region_material->get_options().get_option("x-growth-direction", growth_direction);
		miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
	}
	if (_as->get_options().find_option("y-growth-direction")){
region_material->get_options().get_option("y-growth-direction", growth_direction);
		miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
	}
	if (_as->get_options().find_option("z-growth-direction")){
region_material->get_options().get_option("z-growth-direction", growth_direction);
		miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
		if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
	}  
    

    //-------------------------------------------------------------------------------------------
    
  std::cout << "Material is " << material << std::endl;
      std::cout << "structure is " << structure << std::endl; 
std::cout << "Miller indexes are " << miller  << std::endl;
	
  std::cout << "Ending AtomisticGenerator::do_init() " << std::endl;

}



void 
AtomisticGenerator::set_lattice_type(const std::string lattice_name)
{
#ifdef DEBUG
  std::cerr << "Calling set_lattice_type ... ";
#endif

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
    _prim_vec(3,3) = prim_vec_dir(3,3) * az;

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
  
#ifdef DEBUG
  std::cerr << "done" << endl;
#endif

};


void AtomisticGenerator::set_crystal_basis(const std::string basis_name, const std::string specie1, const std::string specie2, double u)
{
#ifdef DEBUG
  std::cerr << "Calling set_crystal_basis ...";
#endif

  //Set the atomic basis.
  //Vectors are defined in primitive vectors basis.
  _basis_type = basis_name;
  AtomAG tmp;

  assert(_crystal_basis.empty());

  if (_basis_type.compare("cubic") == 0){
    if (specie1.compare("not_specified") == 0) tmp.specie = "Po";
    else tmp.specie = specie1;
     tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("diamond") == 0){
 
    if (specie1.compare("not_specified") == 0) tmp.specie = "C";
    else tmp.specie = specie1;
     tmp.position(1) = -1.0 / 8.0; tmp.position(2) = -1.0/8.0; tmp.position(3) = -1.0/8.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1)=1.0/8.0; tmp.position(2)=1.0/8.0; tmp.position(3)=1.0/8.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("zincblende") == 0){

    if (specie1.compare("not_specified") == 0) tmp.specie = "Ga";
    else tmp.specie = specie1;
    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);
 
    if (specie1.compare("not_specified") == 0) tmp.specie = "As";
    else tmp.specie = specie2;
    tmp.position(1) = 1.0 / 4.0; tmp.position(2) = 1.0 / 4.0; tmp.position(3) = 1.0 / 4.0;
    _crystal_basis.push_back(tmp);
  }

  else if (_basis_type.compare("cristobalite") == 0){

    if (specie1.compare("not_specified") == 0) tmp.specie = "Si";
    else tmp.specie = specie1;
    tmp.position(1) = 1.0 / 8.0; tmp.position(2) = 1.0 / 8.0; tmp.position(3) = 1.0 / 8.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = -1.0 / 8.0; tmp.position(2) = -1.0 / 8.0; tmp.position(3) = -1.0 / 8.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "O";
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

    if (specie1.compare("not_specified") == 0) tmp.specie = "Na";
    else tmp.specie = specie1;
    tmp.position(1) = 0.0; tmp.position(2) = 0.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "Cl";
    else tmp.specie = specie2;
    tmp.position(1) = 1.0 / 2.0; tmp.position(2) =1.0 / 2.0; tmp.position(3) = 1.0 / 2.0;
    _crystal_basis.push_back(tmp);

  }

  else if(_basis_type.compare("wurtzite") == 0){

    if (u == 0.0) u = 3.0 / 8.0;
    if (specie1.compare("not_specified") == 0) tmp.specie = "Ga";
    else tmp.specie = specie1;
    tmp.position(1) = 1.0 / 3.0; tmp.position(2) = 2.0 / 3.0; tmp.position(3) = 0.0;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 2.0 / 3.0; tmp.position(2) = 1.0 / 3.0; tmp.position(3) = 1.0 / 2.0;
    _crystal_basis.push_back(tmp);

    if (specie1.compare("not_specified") == 0) tmp.specie = "N";
    else tmp.specie = specie2;
     tmp.position(1) = 1.0 / 3.0; tmp.position(2) = 2.0 / 3.0; tmp.position(3) = u;
    _crystal_basis.push_back(tmp);

    tmp.position(1) = 2.0 / 3.0; tmp.position(2) = 1.0 / 3.0; tmp.position(3) = 1.0 / 2.0 + u;
    _crystal_basis.push_back(tmp); 
  }

  else if (_basis_type.compare("TiO2") == 0) {

   if (u == 0.0) u = 0.25;
   if (specie1.compare("not_specified") == 0) tmp.specie = "Ti";
   else tmp.specie = specie1;
   tmp.position(1) = -1.0 / 8.0; tmp.position(2) = 5.0 / 8.0; tmp.position(3) = 1.0 / 4.0;
   _crystal_basis.push_back(tmp);

   tmp.position(1) = 1.0 / 8.0; tmp.position(2) = 3.0 / 8.0; tmp.position (3) = 3.0 / 4.0;
   _crystal_basis.push_back(tmp);

   if (specie1.compare("not_specified") == 0) tmp.specie = "O";
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
 std::cerr << "done" << endl;
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






//Some data manipulation function useful only in this class


Tensor2Gen AtomisticGenerator::reciprocal(Tensor2Gen real_basis)
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


int AtomisticGenerator::compare_tol(double a, double b)
{
  //Comparison routine with a tolerance defined as internal constant.
  //If absolute value of difference between a and b is minor than tolerance, 
  //a and b are considered equal
  if (std::fabs(a-b) < tol) return 1;
  else return 0;
};


int AtomisticGenerator::double_to_int_cast_checked(double a)
{
  //Convert a double to the nearest integer, within a certain tolerance
  int n;
  if (std::abs(std::floor(a)-a) < std::abs(std::ceil(a) - a)) n = int(std::floor(a));
  else n = int(std::ceil(a));
  assert (std::abs(double(n) - a) < tol);
  return n;
};


double AtomisticGenerator::double_to_int_value_checked(double a)
{
  //Gives the double number equal to the integer nearest to a, within a certain tolerance
  double b;
  b = double(double_to_int_cast_checked(a));
  return b; 
};


void AtomisticGenerator::double_to_int_value_checked(Tensor1& a)
{
  double tmp;;

  tmp = double_to_int_value_checked(a(1)); a(1) = tmp;
  tmp = double_to_int_value_checked(a(2)); a(2) = tmp;
  tmp = double_to_int_value_checked(a(3)); a(3) = tmp;
};


int AtomisticGenerator::gcd(int a, int b)
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


Tensor1 AtomisticGenerator::reduce_vector(Tensor1 v)
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
  double min_value;
  int i;
  Tensor1 a_tmp;
  Tensor1 a_nonzero;

  if (norm(a) < tol) {a_tmp(1) = 0.0; a_tmp(2) = 0.0; a_tmp(3) = 0.0; a = a_tmp; }

  else {
 
  a_nonzero = a;
  if (compare_tol(a(1),0)) a_nonzero(1) = std::max(std::abs(a(2)),std::abs(a(3)));
  if (compare_tol(a(2),0)) a_nonzero(2) = std::max(std::abs(a(1)),std::abs(a(3)));
  if (compare_tol(a(3),0)) a_nonzero(3) = std::max(std::abs(a(2)),std::abs(a(1)));
 
  min_value = 0;
  //Select minimum value in vector a (0 not allowed from previous constrain)
  min_value = std::min(std::abs(a_nonzero(1)),std::min(std::abs(a_nonzero(2)),std::abs(a_nonzero(3))));

  assert((min_value != 0) && (min_value > 0));

  a_tmp=a * (1 / min_value);
  i = 0;
  do{
    i = i + 1;
    a = a_tmp * i;

  }while ((std::abs(ceil(a(1)) - std::floor(a(1))) > (2 * tol)) && (std::abs(ceil(a(2)) - std::floor(a(2))) > (2 * tol)) && (std::abs(std::ceil(a(3)) - std::floor(a(3))) > (2 * tol)));
  a = reduce_vector(a);
  }

};


void AtomisticGenerator::scale_to_int(Tensor2Gen& a)
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


Tensor2Gen AtomisticGenerator::reduce_vector(Tensor2Gen a)
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

