
#include "DDsemiconductor.h"

using namespace std;
//---------------------------------------------------------------------------------------------//
DDsemiconductor::DDsemiconductor()
{
  strain = Tensor2Sym(0);
  Ev=0; 
  energy_cutoff=1.0; //1eV default value
  strained = false;
}


//--------------------------------------------------------------------------------------------//

DDsemiconductor::DDsemiconductor(const double Ev_1,const  Tensor2Sym& strain_1, const double energy_cutoff_1)
{
  Ev = Ev_1;
  if (norm( strain_1 ) > 1e-5 ) 
    {
      strain = strain_1;
      strained = true;
    }
  else
    {
      strain = Tensor2Sym(0);
      strained = false;
    }
  
  energy_cutoff = energy_cutoff_1;
}
//----------------------------------------------------------------------------------------------//
void DDsemiconductor::set_strain(const Tensor2Sym& strain_1)
{
  
 
  if (norm( strain_1 ) > 1e-5 ) 
    {
      strain = strain_1;
      strained = true;
    }
  else
    {
      strain = Tensor2Sym(0);
      strained = false;
    }
}
//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_conduction_band_energy_mass(void) const
{
  // const std::vector<DDsemiconductor::band_extremum>&  result;

  // result = &conduction_band;

  return(conduction_band);

}
//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_valence_band_energy_mass(void) const
{
  
  return(valence_band);

  
  

}


//---------------------------------------------------------------------------------------------//
DDsemiconductor::~DDsemiconductor (void)
{

}
//---------------------------------------------------------------------------------------------//
