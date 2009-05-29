// $Id$


#include "WzPyroPolarization.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"



WzPyroPolarization*
WzPyroPolarization::create(void)
{
  return new WzPyroPolarization();
}



void
WzPyroPolarization::read_database(void)
{
  Database& db = get_database();
  db.set_section("pyroelectricity");

  _Pz = db.get("Pz", _Pz, true);
}


void
WzPyroPolarization::read_database_alloy(void)
{
  Database& db = get_database();
  db.set_section("pyroelectricity");

  _Pz_bow = db.get("bow_Pz", _Pz_bow);
}



void
WzPyroPolarization::do_init(void)
{
  if (has_parameter("Pz"))
  {
    _Pz = get_option("Pz", _Pz);
    _Pz_bow = 0.0;
  }

  Tensor1 pol(0);
  pol(3) = _Pz;

  set_polarization((get_material()->get_rotated_crystal()).RotMatrix * pol);

  PyroPolarization::do_init();

}


void
WzPyroPolarization::do_init_alloy(const PhysicalModelInterface *comp_A,
        const PhysicalModelInterface *comp_B, double xa)
{
  const WzPyroPolarization* tempA = dynamic_cast<const WzPyroPolarization*>(comp_A);
  const WzPyroPolarization* tempB = dynamic_cast<const WzPyroPolarization*>(comp_B);

  _Pz_bow = get_option("Pz_bow", _Pz_bow);

  _Pz = alloy(tempA->_Pz, tempB->_Pz, xa, _Pz_bow);

  Tensor1 pol(0);
  pol(3) = _Pz;

  set_polarization((get_material()->get_rotated_crystal()).RotMatrix * pol);
}


PhysicalModelInterface*
WzPyroPolarization::create_new(void) const
{
  return new WzPyroPolarization();
}
