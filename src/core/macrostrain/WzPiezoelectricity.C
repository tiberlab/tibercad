// $Id$

#include "WzPiezoelectricity.h"
#include "Database.h"

WzPiezoelectricity::WzPiezoelectricity(const ModelOptions& options) : Piezoelectricity(options)
{
  set_moduli(0,  0,  0);
}

//------------------------------------------------------------//



//-----------------------------------------------------------//

void WzPiezoelectricity::read_database(void)
{

  const Database& db = get_database();
  db.set_section("piezoelectricity");

  e33 = db.get("e33", 0.0);

  e31 = db.get("e31", 0.0);

  e15 = db.get("e15", 0.0);

}


//------------------------------------------------------------//
void WzPiezoelectricity::do_init(void)
{

   get_parameter("e33", e33);

   get_parameter("e31", e31);

   get_parameter("e15", e15);

}


//-----------------------------------------------------------//
void WzPiezoelectricity:: set_moduli(double  e33_i, double e31_i, double e15_i)
{

  e33 = e33_i;
  e31 = e31_i;
  e15 = e15_i;

}

//------------------------------------------------------------//
PhysicalModelInterface* WzPiezoelectricity::create_new(void) const
{
  return new WzPiezoelectricity(get_options());
}



//------------------------------------------------------------//

Tensor1 WzPiezoelectricity::get_polariz_cryst(Tensor2Sym& strain_cryst)
{



  Tensor1 polariz;

  polariz(1) = 2*e15*strain_cryst(3,1);

  polariz(2) = 2*e15*strain_cryst(3,2);

  polariz(3) = e31*strain_cryst(1,1) + e31*strain_cryst(2,2) + e33*strain_cryst(3,3);


  return(polariz);

}

//-----------------------------------------------------------//

void  WzPiezoelectricity::calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const
{
  r = 0;
  r(3,1) = f(1) * e15;
  r(3,2) = f(2) * e15;
  r(1,1) = f(3) * e31;
  r(2,2) = f(3) * e31;
  r(3,3) = f(3) * e33;
}
