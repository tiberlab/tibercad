#include "ZbPiezoelectricity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"

ZbPiezoelectricity::ZbPiezoelectricity() : Piezoelectricity()
{
  e14 = 0;
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

void ZbPiezoelectricity::read_database ()
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());
  e14 = data("e14", 0.0);
}


//---------------------------------------------------------------//

void ZbPiezoelectricity::do_init(void)
{
   ModelOptions & options = get_options ();
   e14 = options.get_option("e14", e14);
  
}

//---------------------------------------------------------------//
PhysicalModelInterface* ZbPiezoelectricity::create_new(void) const
{
  return (new ZbPiezoelectricity());
}

//---------------------------------------------------------------//

void  ZbPiezoelectricity::copy_from (const PhysicalModelInterface* rhs)
{
  const ZbPiezoelectricity* temp = dynamic_cast<const ZbPiezoelectricity*>(rhs);
  e14 = temp->e14;

}

//---------------------------------------------------------------//

void ZbPiezoelectricity::calculate_VCA(const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

   const ZbPiezoelectricity* A = dynamic_cast<const ZbPiezoelectricity*>(comp_A);

   const ZbPiezoelectricity* B = dynamic_cast<const ZbPiezoelectricity*>(comp_B);

   e14 = alloy(A->e14, B->e14, xa);

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
