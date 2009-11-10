// $Id$

#include "ZbStiffness.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"


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

  Database& db = get_database();
  db.set_section("elasticity");

  double c11 = db.get("C11", 0.0, true);
  double c12 = db.get("C12", 0.0, true);
  double c44 = db.get("C44", 0.0, true);

  set_moduli(c11, c12, c44);



}

//-----------------------------------------------//
void ZbStiffness::do_init ( )
{


  double c12 = get_option("C12", C_cr(2,2,1,1));
  double c11 = get_option("C11", C_cr(1,1,1,1));
  double c44 = get_option("C44", C_cr(3,2,3,2));

  set_moduli( c11, c12,  c44);

  Material*   mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calc_system(cr.RotMatrix);

}
