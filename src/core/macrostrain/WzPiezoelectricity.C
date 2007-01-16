#include "WzPiezoelectricity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"

WzPiezoelectricity::WzPiezoelectricity() : Piezoelectricity()
{
  set_moduli(0,  0,  0,  0);
}

//------------------------------------------------------------//

WzPiezoelectricity::WzPiezoelectricity(double  e33, double e31, double e15, double Pz) : Piezoelectricity()
{

  set_moduli(e33,  e31,  e15,  Pz);

}


//-----------------------------------------------------------//

void WzPiezoelectricity::read_database ( )
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  e33 = data("e33", 0.0);
  
  e31 = data("e31", 0.0);

  e15 = data("e15", 0.0);

   Pz = data("Pz", 0.0);
}

//------------------------------------------------------------//
void WzPiezoelectricity::do_init(void)
{
  ModelOptions & options = get_options ();

   e33 = options.get_option("e33",e33);
  
   e31 = options.get_option("e31",e31);

   e15 = options.get_option("e15",e15);

   Pz = options.get_option("Pz", 0.0);

   

}


//-----------------------------------------------------------//
void WzPiezoelectricity:: set_moduli(double  e33_i, double e31_i, double e15_i, double Pz_i)
{

  e33 = e33_i;
  e31 = e31_i;
  e15 = e15_i;
  Pz  = Pz_i;

}

//------------------------------------------------------------//
PhysicalModelInterface* WzPiezoelectricity::create_new(void) const
{
  return new WzPiezoelectricity();
}

//------------------------------------------------------------//

void WzPiezoelectricity::copy_from (const PhysicalModelInterface *rhs)
{
  
  const WzPiezoelectricity* temp = dynamic_cast<const WzPiezoelectricity*> (rhs);
  e33 = temp->e33;
  e31 = temp->e31;
  e15 = temp->e15;
  Pz  = temp->Pz;

} 

//------------------------------------------------------------//

void WzPiezoelectricity:: calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
   const WzPiezoelectricity* tempA = dynamic_cast<const WzPiezoelectricity*> (comp_A);

   const WzPiezoelectricity* tempB = dynamic_cast<const WzPiezoelectricity*> (comp_B);

   e33 = alloy(tempA->e33, tempB->e33, xa);
   e31 = alloy(tempA->e31, tempB->e31, xa);
   e15 = alloy(tempA->e15, tempB->e15, xa);
   Pz  = alloy(tempA->Pz , tempB->Pz , xa);



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
