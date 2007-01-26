using namespace std;
#include "SBZbCondBandBulkHamiltonian.h"
#include "getpot.h"
#include "Alloy.h"
#include "ZbSemiconductor.h"



 
//===========================================================================//
void SBZbCondBandBulkHamiltonian::do_init()
{

  SBCondBandBulkHamiltonian::do_init();

  ZbSemiconductor::ZbDDparameters& par = (dynamic_cast<ZbSemiconductor*> (semiconductor)) -> get_parameters();

  zb_par = &par;

  const ModelOptions& opt =  get_options ();

  min_name = opt.get_option("minimum_name", "Gamma");
  
  if (min_name == "Gamma")
    {
      imass = Tensor2Sym(0);
      imass(1,1) = 1.0/par.m_G;
      imass(2,2) = 1.0/par.m_G;
      imass(3,3) = 1.0/par.m_G;

      edge = ( par.Ev + ((1.0/3.0) * par.delta ) + par.EgGamma)/Hartree;

      kp_bands.resize(1,0);

    }
  
  calculate_Hamiltonian_gen();


}


//===========================================================================//
void SBZbCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
 
  //now strain and potential
  if (min_name == "Gamma")
    {
       Hamiltonian[0][0].constant =  Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
	 + (zb_par->a_c) * trace(strain_crystal )/Hartree;

      
    }
}
