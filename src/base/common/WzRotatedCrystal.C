#include "WzRotatedCrystal.h"

WzRotatedCrystal::WzRotatedCrystal() : RotatedCrystal()
{
  a_lat = 0;
  c_lat = 0;

  set_xyz_mil_direction("x",  1, 0, -1, 0); 
  set_xyz_mil_direction("y", -1, 2, -1, 0);
  set_xyz_mil_direction("z", 0, 0, 0, 1); 



}

//=======================================================//


WzRotatedCrystal::WzRotatedCrystal(const double a, const double c)
{

  set_lat_const(a, c);

}
 
//======================================================//


void WzRotatedCrystal::set_lat_const(const double a, const double c)
{

  a_lat = a;
  c_lat = c;
}

//=======================================================//

void WzRotatedCrystal::calculate_lat_consts()
{

  
  //Bravais vectors 
  Tensor1  Rx;
  Tensor1  Ry;  
  Tensor1  Rz;
  
  // Reciprocal  basis
  Tensor1 rec_basis1;
  Tensor1 rec_basis2;
  Tensor1 rec_basis3;




   Rx = Tensor1(0); Rx(1) = a_lat;
   Ry = Tensor1(0); Ry(1) = -0.5*a_lat; Ry(2) = std::sqrt(3.0)/2.0*a_lat;
   Rz = Tensor1(0); Rz(3) = c_lat;


   const double volume =Rx *vectorProduct(Ry, Rz);

   rec_basis1 = vectorProduct(Ry, Rz)/volume;      
   rec_basis2 = vectorProduct(Rz, Rx)/volume;
   rec_basis3 = vectorProduct(Rx, Ry)/volume;


   lat_const_calc[0] = 1.0/norm((x_miller[0] - x_miller[2]) *rec_basis1  + (x_miller[1]-x_miller[2]) *  rec_basis2
				   + x_miller[3]* rec_basis3);

   lat_const_calc[1] = 1.0/norm((y_miller[0] - y_miller[2]) *rec_basis1  + (y_miller[1]-y_miller[2]) *  rec_basis2
				   + y_miller[3]* rec_basis3);
      
      
   lat_const_calc[2] = 1.0/norm((z_miller[0] - z_miller[2]) *rec_basis1  + (z_miller[1]-z_miller[2]) *  rec_basis2
				   + z_miller[3]* rec_basis3);




}

//======================================================================//

void WzRotatedCrystal::calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil)
{
  assert(vec_x_mil.size()==vec_y_mil.size());
  assert(vec_x_mil.size()==4);

  Tensor1 vec_x;
  Tensor1 vec_y;
  


 // Miller basis
  Tensor1 mil1;
  Tensor1 mil2;
  Tensor1 mil3;


  //assign proncipal directions
      
  mil1=Tensor1(0); mil1(1) = a_lat;
  
  mil2=Tensor1(0); mil2(1) = -0.5*a_lat; mil2(2) = std::sqrt(3.0)/2.0*a_lat;
  
  mil3=Tensor1(0); mil3(3) = c_lat;
  

  //convert from miller indexes to vectors 
      
  vec_x =(x_miller[0] - x_miller[2]) * mil1 + (x_miller[1]-x_miller[2]) * mil2 + x_miller[3]*mil3;
  vec_y =(y_miller[0] - y_miller[2]) * mil1 + (y_miller[1]-y_miller[2]) * mil2 + y_miller[3]*mil3;

  // calculate rotation matrix
  calculate_rot_matrix(vec_x, vec_y);



}

//============================================================================//

void WzRotatedCrystal::set_xyz_mil_direction(std::string dir, int h, int k, int i, int l)
{

  assert( (dir.compare(std::string("x"))==0) || (dir.compare(std::string("y"))==0) || (dir.compare(std::string("z")) ==0) );
  assert((h + k +i == 0));

  if (dir.compare("x")==0) {
    x_miller.clear();

    x_miller.push_back(h);
    x_miller.push_back(k);
    x_miller.push_back(i);
    x_miller.push_back(l);
  };

  if (dir.compare("y")==0) {
    y_miller.clear();


    y_miller.push_back(h);
    y_miller.push_back(k);
    y_miller.push_back(i);
    y_miller.push_back(l);
  };

  
  if (dir.compare("z")==0) {
    z_miller.clear();

    z_miller.push_back(h);
    z_miller.push_back(k);
    z_miller.push_back(i);
    z_miller.push_back(l);
  };

}
//====================================================//
void WzRotatedCrystal::read_database ( )
{


}

//===================================================//

PhysicalModelInterface* WzRotatedCrystal::create_new(void) const
{
  return ( new WzRotatedCrystal() );

}

//===================================================//

void WzRotatedCrystal::copy_from (const PhysicalModelInterface *rhs) 
{

  const WzRotatedCrystal* temp = dynamic_cast<const WzRotatedCrystal*> (rhs);


  a_lat = temp->a_lat;

  c_lat = temp->c_lat;


  x_miller = temp->x_miller;

  y_miller = temp->y_miller;

  z_miller = temp->z_miller;

  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);

}

//===================================================//


void WzRotatedCrystal::do_init(void)
{
  ModelOptions & options = get_options ();

  a_lat =  options.get_option("a", a_lat);

  assert(a_lat > 0);

  c_lat =  options.get_option("c", c_lat);

  assert(c_lat > 0);

  options.get_option("x-growth-direction", x_miller);
  options.get_option("y-growth-direction", y_miller);
  options.get_option("z-growth-direction", z_miller);
  

  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);


}

//===================================================//

void WzRotatedCrystal::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

 

  const WzRotatedCrystal* modA = dynamic_cast<const WzRotatedCrystal*> (comp_A);

  const WzRotatedCrystal* modB = dynamic_cast<const WzRotatedCrystal*> (comp_B);


  x_miller = modA->x_miller;

  y_miller = modA->y_miller;

  z_miller = modA->z_miller;


  a_lat = alloy(modA->a_lat, modB->a_lat, xa);

  
  c_lat = alloy(modA->c_lat, modB->c_lat, xa);


  


  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);
  

}
