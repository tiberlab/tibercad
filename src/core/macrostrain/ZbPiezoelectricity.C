// $Id$

#include "ZbPiezoelectricity.h"
#include "Material.h"
#include "Database.h"

ZbPiezoelectricity::ZbPiezoelectricity(const ModelOptions& options) : Piezoelectricity(options)
{
  e14 = 0;
}



//---------------------------------------------------------------//


void ZbPiezoelectricity::set_piezo_module(double e)
{
  e14 = e;
}

//---------------------------------------------------------------//

void ZbPiezoelectricity::read_database ()
{

  Database& db = get_database();
  db.set_section("piezoelectricity");

  e14 = db.get("e14", 0.0);
}



//---------------------------------------------------------------//

void ZbPiezoelectricity::do_init(void)
{

   get_parameter("e14", e14);

}

//---------------------------------------------------------------//
PhysicalModelInterface* ZbPiezoelectricity::create_new(void) const
{
  return (new ZbPiezoelectricity(get_options()));
}



//---------------------------------------------------------------//

Tensor1  ZbPiezoelectricity::get_polariz_cryst(Tensor2Sym& strain_cryst)
{

  Tensor1 polariz;

  polariz(1) = 2*e14*strain_cryst(3,2);
  polariz(2) = 2*e14*strain_cryst(3,1);
  polariz(3) = 2*e14*strain_cryst(2,1);

  return(polariz);
}
//-----------------------------------------------------------------//
void  ZbPiezoelectricity::calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const
{
  r = 0;
  r(3,2) = f(1) * e14;
  r(3,1) = f(2) * e14;
  r(2,1) = f(3) * e14;
}

