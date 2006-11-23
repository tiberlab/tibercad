#include "RotatedCrystal.h"

RotatedCrystal::RotatedCrystal() : PhysicalProperties("Rotated Crystal")
{
  
}

//===============================================================//


void RotatedCrystal::calculate_rot_matrix(const Tensor1& vec_x, const Tensor1& vec_y)
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


//===================================================================//

Tensor2Sym RotatedCrystal::get_const_eps0(double lat_cont_substrate[3], Tensor2Sym& eps0_var_log)
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

//========================================================================//

Tensor2Sym RotatedCrystal::get_var_eps0( std::string  var_name)
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
 
//=======================================================================//

void RotatedCrystal::get_lat_const(double lat_const[3])
{
  for (int i = 0; i <=2; i++) lat_const[i] = lat_const_calc[i];
}


//=======================================================================//


Tensor2Sym RotatedCrystal::get_eps0(double lat_cont_substrate[3])
{
  Tensor2Sym eps0;
  eps0 = Tensor2Sym(0);
  for (int i = 1; i <=3; i++)
    {
      eps0(i,i) = (lat_cont_substrate[i-1] - lat_const_calc[i-1])/lat_const_calc[i-1];
    }

  return(eps0); 

}


//========================================================================//
