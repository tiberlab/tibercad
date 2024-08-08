// $Id$

#include "RotatedCrystal.h"
#include "ModelErrorException.h"
#include "Material.h"


RotatedCrystal::RotatedCrystal(const ModelOptions& options)
 : PhysicalModel(options)
{
  x_miller.resize(3, 0);
  x_miller[0] = 1;
  y_miller.resize(3, 0);
  y_miller[1] = 1;
  z_miller.resize(3, 0);
  z_miller[2] = 1;
}



RotatedCrystal*
RotatedCrystal::create(const Material* owner,
    const ModelOptions& options)
{
  std::string structure = owner->get_structure();

  RotatedCrystal* rc = NULL;

  if (structure != "zb" && structure != "wz")
    structure = "gen";
    
  if ((structure == "am") || (structure == "amorphous"))
    structure = "gen";
  
  rc = PhysicalModel::create<RotatedCrystal>(
      "cryst_" + structure, owner, options);

  if (rc == NULL)
  {
    std::string msg("No such crystal structure known: ");
    msg += structure;
    throw ModelErrorException(msg);
  }

  rc->set_material(owner);

  return rc;
}






//===============================================================//


void RotatedCrystal::calculate_rot_matrix(const Tensor1& vec_x, const Tensor1& vec_y)
{
 
  //calculate 3rd vector in the basis
  const Tensor1 vec_z = vectorProduct(vec_x, vec_y);

  double scal = fabs(vec_x * vec_y);
  if (scal > 1e-6)
  {
    throw ModelErrorException("Provided x and y directions are non-orthogonal!");
  }

  assert(norm(vec_x) > 1e-10);
  assert(norm(vec_y) > 1e-10);
  assert(norm(vec_z) > 1e-10);

  //create rotation matrix


  for ( int i = 1; i <=3; i++ )  RotMatrix(1,i) = vec_x(i)/norm(vec_x);
  for ( int i = 1; i <=3; i++ )  RotMatrix(2,i) = vec_y(i)/norm(vec_y);
  for ( int i = 1; i <=3; i++ )  RotMatrix(3,i) = vec_z(i)/norm(vec_z);

}


//===================================================================//

Tensor2Sym RotatedCrystal::get_const_eps0(double lat_cont_substrate[3], Tensor2Sym& eps0_var_log) const
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

Tensor2Sym RotatedCrystal::get_var_eps0( std::string  var_name) const
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

void RotatedCrystal::get_lat_const(double lat_const[3]) const
{
  for (int i = 0; i <=2; i++) lat_const[i] = lat_const_calc[i];
}


//=======================================================================//


Tensor2Sym RotatedCrystal::get_eps0(double lat_cont_substrate[3]) const
{
  Tensor2Sym eps0;
  eps0 = Tensor2Sym(0);
  for (int i = 1; i <=3; i++)
    {
      eps0(i,i) = (lat_cont_substrate[i-1] - lat_const_calc[i-1])/lat_const_calc[i-1];
    }

  return(eps0); 

}


void RotatedCrystal::calculate_euler_angles(void)
{
  // calculate Euler angles

  // formulas taken from en.wikipedia.org/wiki/Euler_angles
  // and https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf

  // Note: they are calculated from the transpose of RotMatrix
  _beta = acos(RotMatrix(3, 3));

  if (abs(RotMatrix(3,3)) < (1.0 - 1e-6))
  {

    double sb1 = sin(_beta);
    _alpha = atan2(-RotMatrix(2,3)/sb1, RotMatrix(1,3)/sb1);


    if (abs(_alpha) > M_PI_2)
    {
      double sb2 = -sb1;
      _alpha = atan2(-RotMatrix(2,3)/sb2, RotMatrix(1,3)/sb2);
      _gamma = atan2(RotMatrix(3,2)/sb2, RotMatrix(3,1)/sb2);
      _beta = -_beta;
    }
    else
      _gamma = atan2(RotMatrix(3,2)/sb1, RotMatrix(3,1)/sb1);
  }
  else // R33 = +/-1
  {
    // set gamma = 0 arbitrarily
    _gamma = 0;

    if (RotMatrix(3,3) < 0) // R33 = -1
    {
      _alpha = atan2(RotMatrix(2,1), RotMatrix(2,2));
    }
    else
    {
      _alpha = atan2(RotMatrix(1,2), RotMatrix(1,1));
    }


  }
}



void RotatedCrystal::do_init(void)
{
  Tensor1 vec_x;
  Tensor1 vec_y;

  vec_x(1) = x_miller[0];
  vec_x(2) = x_miller[1];
  vec_x(3) = x_miller[2];

  vec_y(1) = y_miller[0];
  vec_y(2) = y_miller[1];
  vec_y(3) = y_miller[2];

  calculate_rot_matrix(vec_x, vec_y);

  calculate_euler_angles();


}

//========================================================================//
