#include "ZbStiffness.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"

//---------------------------------------------------//

ZbStiffness::ZbStiffness() :Stiffness()
{
set_moduli(0, 0,  0);
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
void ZbStiffness::read_database ( )
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  double c11 = data("C11", 0);
  double c12 = data("C12", 0);
  double c44 = data("C44", 0);

  set_moduli( c11, c12, c44);

  
  
}

//-----------------------------------------------//
void ZbStiffness::do_init ( )
{
  ModelOptions & options = get_options ();
  
  double c12 = options.get_option ("C12", C_cr(2,2,1,1)); 
  double c11 = options.get_option ("C11", C_cr(1,1,1,1));
  double c44 = options.get_option ("C44", C_cr(3,2,3,2));

  set_moduli( c11, c12,  c44);

}
