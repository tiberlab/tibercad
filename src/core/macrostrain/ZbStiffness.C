#include "ZbStiffness.h"

//---------------------------------------------------//

ZbStiffness::ZbStiffness() :Stiffness()
{

}
//---------------------------------------------------//

ZbStiffness::ZbStiffness(double c11, double c12, double c44) : Stiffness()
{   

  set_moduli(c11,  c12,  c44);

} 
 
//----------------------------------------------------//

void ZbStiffness::set_moduli(double c11, double c12, double c44)
{

  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c11;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c12;
  C_cr(3,3,2,2) = c12;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;

  C_calc = C_cr;

}

//----------------------------------------------//
void ZbStiffness::read_database (const Dummy &db)
{
 
}
