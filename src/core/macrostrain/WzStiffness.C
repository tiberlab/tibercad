// $Id$

#include "WzStiffness.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"

//---------------------------------------------------//

WzStiffness::WzStiffness(const ModelOptions& options) :Stiffness(options)
{
  set_moduli(0,0,0,0,0);
}
//---------------------------------------------------//


void WzStiffness::set_moduli(double c11, double c12, double c13, double c33, double c44)

{
  C_cr(1,1,1,1) = c11;
  C_cr(2,2,2,2) = c11;
  C_cr(3,3,3,3) = c33;
  C_cr(2,2,1,1) = c12;
  C_cr(3,3,1,1) = c13;
  C_cr(3,3,2,2) = c13;
  C_cr(2,1,2,1) = 0.5*(c11-c12);
  C_cr(3,1,3,1) = c44;
  C_cr(3,2,3,2) = c44;

  C_calc = C_cr;

}

//----------------------------------------------//
void WzStiffness::read_database ( )
{

  const Database& db = get_database();
  db.set_section("elasticity");

  double c11 = db.get("C11", 0.0, true);
  double c12 = db.get("C12", 0.0, true);
  double c13 = db.get("C13", 0.0, true);
  double c33 = db.get("C33", 0.0, true);
  double c44 = db.get("C44", 0.0, true);


  set_moduli(c11,  c12,  c13,  c33,  c44);


}


//----------------------------------------------//
void WzStiffness::do_init(void)
{



  double c11 = get_option("C11", C_cr(1,1,1,1));
  double c12 = get_option("C12", C_cr(2,2,1,1));
  double c13 = get_option("C13", C_cr(3,3,1,1));
  double c33 = get_option("C33", C_cr(3,3,3,3));
  double c44 = get_option("C44", C_cr(3,2,3,2));

  set_moduli(c11,  c12,  c13,  c33,  c44);

  const Material*   mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calc_system(cr.RotMatrix);




}



