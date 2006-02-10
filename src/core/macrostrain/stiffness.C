#include <tensor.h>
#include <cmath>
#include <xtensor.h>
using std::cout;

#include "stiffness.h"


//-------------------------------------------------------------------------------------
stiffness::stiffness()
{
  C_cr = 0;
}

//--------------------------------------------------------------------------------------------------
stiffness::stiffness(double c11, double c12, double c44)
{
  // assembles stiffness tensor in crystal system for a zinc-blende crystal
  // c11 = C_{xxxx}, c12=C_{xxyy}, c44=C_{xyxy} (Voigt notation) 
  C_cr = 0;
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c11;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c12;
  C_cr(3,3,2,2) = c12;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;
 }

//--------------------------------------------------------------------------------------------------
stiffness::stiffness(double c11, double c12, double c13, double c33, double c44)
{
  //assembles stiffness tensor in crystal system for a wurtzite crystal
  //z axis is parallel to [0001] direction
  //(x,y,z) - othogonormal system 
  // c11 = C_{xxxx}, c12=C_{xxyy}, c13 = C_{xxzz}, c33 = C_{zzzz}, c44 = C_{yzxz}, (Voigt notation) 
  
  C_cr = 0;
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c33;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c13;
  C_cr(3,3,2,2) = c13;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;
 
 

}

void stiffness::set_moduli(double c11, double c12, double c44)
{
  // assembles stiffness tensor in crystal system for a zinc-blende crystal
  // c11 = C_{xxxx}, c12=C_{xxyy}, c44=C_{xyxy} (Voigt notation) 
  C_cr = 0;
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c11;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c12;
  C_cr(3,3,2,2) = c12;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;
 }

//--------------------------------------------------------------------------------------------------
void stiffness::set_moduli(double c11, double c12, double c13, double c33, double c44)
{
  //assembles stiffness tensor in crystal system for a wurtzite crystal
  //z axis is parallel to [0001] direction
  //(x,y,z) - othogonormal system 
  // c11 = C_{xxxx}, c12=C_{xxyy}, c13 = C_{xxzz}, c33 = C_{zzzz}, c44 = C_{yzxz}, (Voigt notation) 
  
  C_cr = 0;
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c33;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c13;
  C_cr(3,3,2,2) = c13;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;


}
  //-------------------------------------------------------------------------------------
void stiffness::print_crystall_system()
{
  //prints NOT in Voigt notation
  cout << std::setw(12) <<C_cr << "\n";
}
//-------------------------------------------------------------------------------------

void stiffness::rotate_to_calc_system(const Tensor2Gen& RotMatrix)
{
  // generates stiffness matrix in calculation system
 
  C_calc = push_forward(C_cr, RotMatrix);
}


