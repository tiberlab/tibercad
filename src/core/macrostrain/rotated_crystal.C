#include "rotated_crystal.h"

using std::cout;
using std::string;


//-------------------------------------------------------

rotated_crystal::rotated_crystal()
{
  RotMatrix = Tensor2Gen(1) ;
}


//--------------------------------------------------------

rotated_crystal::rotated_crystal(string type, double a)
{
  
  assert(type.compare("cub")==0);
  crystal_type = type;
  RotMatrix = Tensor2Gen(1) ;
  a_lat = a;
  
}

//---------------------------------------------------------

rotated_crystal::rotated_crystal(string type, double a, double c)
{
  assert(type.compare("hex")==0);
  crystal_type = type;
  RotMatrix = Tensor2Gen(1) ;
  a_lat = a;
  c_lat = c;
}

//--------------------------------------------------------

void rotated_crystal:: set_cryst_type(const string type)
{
  const string str1 = "cub";
  const string str2 = "hex";
  assert ((type.compare(str1)==0) || (type.compare(str2)==0)) ;
  crystal_type = type;
  
} 

//-----------------------------------------------------------

string rotated_crystal::get_cryst_type()
{
  return(crystal_type);
}

//-----------------------------------------------------------

void rotated_crystal::set_lat_const(const double a)
{
  //cubic lattice - all constants are the same!
  a_lat = a;
  
}

//------------------------------------------------------------

void rotated_crystal::set_lat_const(const double a, const double c )
{
  //hexagonal lattice - firsttwo constants are equal!
  a_lat = a;
  c_lat = c;
  
}

//-------------------------------------------------------------

void rotated_crystal::set_xyz_mil_direction(string dir, int h, int k, int l)
{
  assert( (dir.compare(string("x"))==0) || (dir.compare(string("y"))==0) || (dir.compare(string("z")) ==0) );

    if (dir.compare("x")==0) {
     
      x_miller.push_back(h);
      x_miller.push_back(k);
      x_miller.push_back(l);
  };

    if (dir.compare("y")==0) {
      y_miller.push_back(h);
      y_miller.push_back(k);
      y_miller.push_back(l);
      
  };

  
    if (dir.compare("z")==0) {
      z_miller.push_back(h);
      z_miller.push_back(k);
      z_miller.push_back(l);
     
  };
}
//-------------------------------------------------------  
void rotated_crystal::set_xyz_mil_direction(string dir, int h, int k, int i, int l)
{ 
  assert( (dir.compare(string("x"))==0) || (dir.compare(string("y"))==0) || (dir.compare(string("z")) ==0) );
  assert((h + k +i == 0));

  if (dir.compare("x")==0) {
    x_miller.push_back(h);
    x_miller.push_back(k);
    x_miller.push_back(i);
    x_miller.push_back(l);
  };

  if (dir.compare("y")==0) {
    y_miller.push_back(h);
    y_miller.push_back(k);
    y_miller.push_back(i);
    y_miller.push_back(l);
  };

  
  if (dir.compare("z")==0) {
    z_miller.push_back(h);
    z_miller.push_back(k);
    z_miller.push_back(i);
    z_miller.push_back(l);
  };

}
//---------------------------------------------------------

void  rotated_crystal::calculate_lat_consts()
{

 


 
  //---------------------------------------------------------
  //Bravais vectors 
  Tensor1  Rx;
  Tensor1  Ry;  
  Tensor1  Rz;
  //---------------------------------------------------------
  // Miller basis
  Tensor1 mil1;
  Tensor1 mil2;
  Tensor1 mil3;

  
  Tensor1 rec_basis1;
  Tensor1 rec_basis2;
  Tensor1 rec_basis3;


  if (crystal_type.compare("cub")==0)
    {
     
      //calculate lattice constants as in-plane distances
  


      lat_const_calc[0] = a_lat/sqrt(double( x_miller[0]*x_miller[0] + x_miller[1]*x_miller[1] + x_miller[2]*x_miller[2]));
      lat_const_calc[1] = a_lat/sqrt(double( y_miller[0]*y_miller[0] + y_miller[1]*y_miller[1] + y_miller[2]*y_miller[2]));
      lat_const_calc[2] = a_lat/sqrt(double( z_miller[0]*z_miller[0] + z_miller[1]*z_miller[1] + z_miller[2]*z_miller[2]));


    }

  if (crystal_type.compare("hex")==0)
    {


      for (short i = 0; i<3; i++)
	lat_const_calc[i] = a_lat;

      if (z_miller[3] != 0) lat_const_calc[2] = c_lat;
      if (y_miller[3] != 0) lat_const_calc[1] = c_lat;
      if (x_miller[3] != 0) lat_const_calc[0] = c_lat;
    
    /*


      mil1=Tensor1(0); mil1(1) = a_lat;
      
      mil2=Tensor1(0); mil2(1) = -0.5*a_lat; mil2(2) = std::sqrt(3.0)/2.0*a_lat;
      
      mil3=Tensor1(0); mil3(3) = c_lat;

      

      Rx =(x_miller[0] - x_miller[2]) * mil1 + (x_miller[1]-x_miller[2]) * mil2 + x_miller[3]*mil3;
      Ry =(y_miller[0] - y_miller[2]) * mil1 + (y_miller[1]-y_miller[2]) * mil2 + y_miller[3]*mil3;
      Rz =(z_miller[0] - z_miller[2]) * mil1 + (z_miller[1]-z_miller[2]) * mil2 + z_miller[3]*mil3;

      
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
       


      */

    }

}


//---------------------------------------------------------------------------------------
Tensor2Sym rotated_crystal::get_eps0(double lat_cont_substrate[3])
{
  Tensor2Sym eps0;
  eps0 = Tensor2Sym(0);
  for (int i = 1; i <=3; i++)
    {
      eps0(i,i) = (lat_cont_substrate[i-1] - lat_const_calc[i-1])/lat_const_calc[i-1];
    }

  return(eps0); 

}
//-------------------------------------------------------------------------------
Tensor2Sym rotated_crystal::get_const_eps0(double lat_cont_substrate[3], Tensor2Sym& eps0_var_log)
{
  Tensor2Sym const_eps0;
  const_eps0 = Tensor2Sym(0);

 


  for (int i = 1; i <=3; i++)
    {
      if (eps0_var_log(i,i) == 1)
	{
	  const_eps0(i,i) = -1.0;
	}
      else
	{
	  
	  const_eps0(i,i) = (lat_cont_substrate[i-1] - lat_const_calc[i-1])/lat_const_calc[i-1];
	}
    }

  return(const_eps0);

}
//--------------------------------------------------------------------------------
Tensor2Sym rotated_crystal::get_var_eps0( std::string  var_name)
{
  Tensor2Sym eps0;
  eps0 = Tensor2Sym(0);
  if (var_name == "ax")
    {    
      eps0(1,1) = 1.0/lat_const_calc[1-1];
    }

  
  if (var_name == "ay")
    {      
      eps0(2,2) = 1.0/lat_const_calc[2-1];
    }

  
  if (var_name == "az")
    {
     
      eps0(3,3) = 1.0/lat_const_calc[3-1];
    }
  
  
  if (var_name == "eps_xy")
    {
      eps0(2,1) = 1;
    }

  if (var_name == "eps_xz")
    {
      eps0(3,1) = 1;
    }

  if (var_name == "eps_yz")
    {
      eps0(3,2) = 1;
    }

  return(eps0);

}



//---------------------------------------------------------------------------------
void rotated_crystal::get_lat_const(double lat_const[3])
{
  for (int i = 0; i <=2; i++) lat_const[i] = lat_const_calc[i];
}

//------------------------------------------------------------

 void rotated_crystal::calculate_rot_matrix(const Tensor1& vec_x, const Tensor1& vec_y)
{
 
  //calculate 3rd vector in the basis
  const Tensor1 vec_z = vectorProduct(vec_x, vec_y) ;

  assert(norm(vec_x) > 1e-10);
  assert(norm(vec_y) > 1e-10);
  assert(norm(vec_z) > 1e-10);

  //create rotation matrix
  for ( int i = 1; i <=3; i++ )  RotMatrix(1,i) = vec_x(i)/norm(vec_x);
  for ( int i = 1; i <=3; i++ )  RotMatrix(2,i) = vec_y(i)/norm(vec_y);
  for ( int i = 1; i <=3; i++ )  RotMatrix(3,i) = vec_z(i)/norm(vec_z);


}

//------------------------------------------------------------------

void rotated_crystal::calculate_rot_matrix(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil)
{
  assert(vec_x_mil.size()==vec_y_mil.size());
  assert(vec_x_mil.size()==3 || vec_x_mil.size()==4);
  Tensor1 vec_x;
  Tensor1 vec_y;
  if (vec_x_mil.size()==3)
    {
      assert(crystal_type.compare( "cub") == 0);
     
      //convert from miller indexes to vectors
      
      vec_x(1) = a_lat* vec_x_mil[0];
      vec_x(2) = a_lat* vec_x_mil[1];
      vec_x(3) = a_lat* vec_x_mil[2];
      
      vec_y(1) = a_lat* vec_y_mil[0];
      vec_y(2) = a_lat* vec_y_mil[1];
      vec_y(3) = a_lat* vec_y_mil[2];
      
      
    }

  if (vec_x_mil.size()==4)
    {
       // Miller basis
      Tensor1 mil1;
      Tensor1 mil2;
      Tensor1 mil3;

      assert(crystal_type.compare( "hex") == 0);
    
      //assign proncipal directions
      
      mil1=Tensor1(0); mil1(1) = a_lat;
  
      mil2=Tensor1(0); mil2(1) = -0.5*a_lat; mil2(2) = std::sqrt(3.0)/2.0*a_lat;
  
      mil3=Tensor1(0); mil3(3) = c_lat;
      
      //convert from miller indexes to vectors 
      
      vec_x =(x_miller[0] - x_miller[2]) * mil1 + (x_miller[1]-x_miller[2]) * mil2 + x_miller[3]*mil3;
      vec_y =(y_miller[0] - y_miller[2]) * mil1 + (y_miller[1]-y_miller[2]) * mil2 + y_miller[3]*mil3;    
    }

  // calculate rotation matrix
      
  calculate_rot_matrix(vec_x, vec_y);

}

//------------------------------------------------------------------


