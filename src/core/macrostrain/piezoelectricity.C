# include "piezoelectricity.h"
//-----------------------------------------------------------------//
Piezoelectricity :: Piezoelectricity()
{
  pyro_const = 0.0;
}


//---------------------------------------------------------------//

void Piezoelectricity :: set_moduli(double  e14)
{
  moduli.clear();
  type = "zb";
  moduli.push_back(e14);
}

//----------------------------------------------------------------//

void Piezoelectricity :: set_moduli(double  e31, double e15, double e33)
{
  moduli.clear();
  type = "wz";
  moduli.push_back(e31);
  moduli.push_back(e15);
  moduli.push_back(e33);
}
//-----------------------------------------------------------------//
void Piezoelectricity :: set_pyro_module(double  p)
{
  pyro_const = p; 
}

//----------------------------------------------------------------//

Tensor1 Piezoelectricity :: get_polariz_cryst(Tensor2Sym& strain_cryst)
{
 
  Tensor1 polariz;

  if (type=="zb")
    {
      double e14 = moduli[0];

     
      polariz(1) = 2*e14*strain_cryst(3,2);
      polariz(2) = 2*e14*strain_cryst(3,1);
      polariz(3) = 2*e14*strain_cryst(2,1);
    }
  //----------------//
  if (type=="wz")
    {
      double e31 = moduli[0];
      double e15 = moduli[1];
      double e33 = moduli[2];

      polariz(1) = 2*e15*strain_cryst(3,1);
      polariz(2) = 2*e15*strain_cryst(2,1);
      polariz(3) = e31*strain_cryst(1,1) + e31*strain_cryst(2,2) + e33*strain_cryst(3,3);

      polariz(3) += pyro_const;

    }

  return(polariz);

}
