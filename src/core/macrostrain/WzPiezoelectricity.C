#include "WzPiezoelectricity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"

WzPiezoelectricity::WzPiezoelectricity() : Piezoelectricity()
{
  set_moduli(0,  0,  0,  0);

  e33_bow = 0;
  e31_bow = 0;
  e15_bow = 0;
  Pz_bow = 0;
}

//------------------------------------------------------------//

WzPiezoelectricity::WzPiezoelectricity(double  e33, double e31, double e15, double Pz) : Piezoelectricity()
{

  set_moduli(e33,  e31,  e15,  Pz);


  e33_bow = 0;
  e31_bow = 0;
  e15_bow = 0;
  Pz_bow = 0;

}


//-----------------------------------------------------------//

void WzPiezoelectricity::read_database(void)
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  e33 = data("e33", 0.0);
  
  e31 = data("e31", 0.0);

  e15 = data("e15", 0.0);

  Pz = data("Pz", 0.0);

}


void
WzPiezoelectricity::read_database_alloy(void)
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  e33_bow = data("bow_e33", 0.0);
  
  e31_bow = data("bow_e31", 0.0);

  e15_bow = data("bow_e15", 0.0);

  Pz_bow = data("bow_Pz", 0.0);

}


//------------------------------------------------------------//
void WzPiezoelectricity::do_init(void)
{

   e33 = get_parameter("e33",e33);
  
   e31 = get_parameter("e31",e31);

   e15 = get_parameter("e15",e15);

   Pz = get_parameter("Pz", Pz);
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

void WzPiezoelectricity:: do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
   const WzPiezoelectricity* tempA = dynamic_cast<const WzPiezoelectricity*> (comp_A);
   const WzPiezoelectricity* tempB = dynamic_cast<const WzPiezoelectricity*> (comp_B);


   e33_bow = get_parameter("bow_e33",e33);
   e31_bow = get_parameter("bow_e31",e31);
   e15_bow = get_parameter("bow_e15",e15);
   Pz_bow = get_parameter("bow_Pz", Pz);


   e33 = alloy(tempA->e33, tempB->e33, xa, e33_bow );
   e31 = alloy(tempA->e31, tempB->e31, xa, e31_bow );
   e15 = alloy(tempA->e15, tempB->e15, xa, e15_bow );
   Pz  = alloy(tempA->Pz , tempB->Pz , xa, Pz_bow  );

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
