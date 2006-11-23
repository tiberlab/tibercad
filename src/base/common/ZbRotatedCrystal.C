#include "ZbRotatedCrystal.h"

ZbRotatedCrystal::ZbRotatedCrystal() : RotatedCrystal()
{

}

//===============================================================//


ZbRotatedCrystal::ZbRotatedCrystal(const double a) : RotatedCrystal()
{
  set_lat_const(a);
}

//===============================================================//


void ZbRotatedCrystal::set_lat_const(const double a) 
{
  a_lat = a;
}

//===============================================================//

void ZbRotatedCrystal::set_xyz_mil_direction(std::string dir, int h, int k, int l)
{
  assert( (dir.compare(std::string("x"))==0) || (dir.compare(std::string("y"))==0) || (dir.compare(std::string("z")) ==0) );

   if (dir.compare("x")==0) {
     x_miller.clear();
 
     x_miller.push_back(h);
     x_miller.push_back(k);
     x_miller.push_back(l);
  };

    if (dir.compare("y")==0) {
      y_miller.clear();

      y_miller.push_back(h);
      y_miller.push_back(k);
      y_miller.push_back(l);
      
  };

  
    if (dir.compare("z")==0) {
      z_miller.clear();

      z_miller.push_back(h);
      z_miller.push_back(k);
      z_miller.push_back(l);
     
  };

}

//===============================================================//

void ZbRotatedCrystal::calculate_lat_consts()
{

  

 lat_const_calc[0] = a_lat/sqrt(double( x_miller[0]*x_miller[0] + x_miller[1]*x_miller[1] + x_miller[2]*x_miller[2]));
 lat_const_calc[1] = a_lat/sqrt(double( y_miller[0]*y_miller[0] + y_miller[1]*y_miller[1] + y_miller[2]*y_miller[2]));
 lat_const_calc[2] = a_lat/sqrt(double( z_miller[0]*z_miller[0] + z_miller[1]*z_miller[1] + z_miller[2]*z_miller[2]));


 std::cerr << "lat_const_calc  " << lat_const_calc[0] << "   " << lat_const_calc[1] << "    " <<  lat_const_calc[2] << "\n";
 

}

//===============================================================//

void ZbRotatedCrystal::calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil)
{

 assert(vec_x_mil.size()==vec_y_mil.size());
 assert(vec_x_mil.size()==3);


 Tensor1 vec_x;
 Tensor1 vec_y;

 //convert from miller indexes to vectors

 vec_x(1) = a_lat* vec_x_mil[0];
 vec_x(2) = a_lat* vec_x_mil[1];
 vec_x(3) = a_lat* vec_x_mil[2];
      
 vec_y(1) = a_lat* vec_y_mil[0];
 vec_y(2) = a_lat* vec_y_mil[1];
 vec_y(3) = a_lat* vec_y_mil[2];


 calculate_rot_matrix(vec_x, vec_y);



}


//====================================================//
void ZbRotatedCrystal::read_database (Dummy &db)
{

   
}
