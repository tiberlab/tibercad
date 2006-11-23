#include "WzPiezoelectricity.h"

WzPiezoelectricity::WzPiezoelectricity() : Piezoelectricity()
{

}

//------------------------------------------------------------//

WzPiezoelectricity::WzPiezoelectricity(double  e33, double e31, double e15, double Pz) : Piezoelectricity()
{

  set_moduli(e33,  e31,  e15,  Pz);

}


//-----------------------------------------------------------//

void WzPiezoelectricity::read_database (const Dummy &db)
{


}

//------------------------------------------------------------//

void WzPiezoelectricity:: set_moduli(double  e33_i, double e31_i, double e15_i, double Pz_i)
{

  e33 = e33_i;
  e31 = e31_i;
  e15 = e15_i;
  Pz  = Pz_i;

}

//------------------------------------------------------------//

Tensor1 WzPiezoelectricity::get_polariz_cryst(Tensor2Sym& strain_cryst)
{

  Tensor1 polariz;

  polariz(1) = 2*e15*strain_cryst(3,1);
  polariz(2) = 2*e15*strain_cryst(2,1);
  polariz(3) = e31*strain_cryst(1,1) + e31*strain_cryst(2,2) + e33*strain_cryst(3,3);

  polariz(3) += Pz;

}
