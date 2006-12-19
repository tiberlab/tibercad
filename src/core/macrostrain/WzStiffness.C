#include "WzStiffness.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"

//---------------------------------------------------//

WzStiffness::WzStiffness() :Stiffness()
{
  set_moduli(0,0,0,0,0);
}  
//---------------------------------------------------//

WzStiffness::WzStiffness(double c11, double c12, double c13, double c33, double c44) :Stiffness()
{
  set_moduli(c11, c12,  c13,  c33,  c44);
}

//-----------------------------------------------//

void WzStiffness::set_moduli(double c11, double c12, double c13, double c33, double c44)

{
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c33;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c13;
  C_cr(3,3,2,2) = c13;
  C_cr(2,1,2,1) = c44;
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;

  C_calc = C_cr;

}

//----------------------------------------------//
void WzStiffness::read_database ( )
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());
  double c11 = data ("C11", 0);
  double c12 = data ("C12", 0);
  double c13 = data ("C13", 0);
  double c33 = data ("C33", 0);
  double c44 = data ("C44", 0);

  set_moduli(c11,  c12,  c13,  c33,  c44);


}


//----------------------------------------------//
void WzStiffness::do_init(void)
{
   ModelOptions & options = get_options ();
   double c11 = options.get_option ("C11", C_cr(1,1,1,1));
   double c12 = options.get_option ("C12", C_cr(2,2,1,1));
   double c13 = options.get_option ("C13", C_cr(3,3,1,1));
   double c33 = options.get_option ("C33", C_cr(3,3,3,3));
   double c44 = options.get_option ("C44", C_cr(3,2,3,2));

   set_moduli(c11,  c12,  c13,  c33,  c44);

}



