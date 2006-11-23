#include "ZbPiezoelectricity.h"

ZbPiezoelectricity::ZbPiezoelectricity() : Piezoelectricity()
{

}

//---------------------------------------------------------------//

ZbPiezoelectricity::ZbPiezoelectricity(double e14) : Piezoelectricity()
{
  set_piezo_module(e14);
}

//---------------------------------------------------------------//


void ZbPiezoelectricity::set_piezo_module(double e)
{
  e14 = e;
}

//---------------------------------------------------------------//

void ZbPiezoelectricity::read_database (const Dummy &db)
{

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
