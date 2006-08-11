using namespace std;
#include "SBZbCondBandBulkHamiltonian.h"
#include "getpot.h"
#include "Alloy.h"


SBZbCondBandBulkHamiltonian::SBZbCondBandBulkHamiltonian( ):SBbulkHamiltonian()
{
  kp_bands.resize(1,0);
  kp_bands_map.insert(make_pair(0,0));
  min_name = "Gamma";
}

//=======================================================================//
void  SBZbCondBandBulkHamiltonian::calculate_edge_and_mass()
{
  //--------------------------------------------------
  //if there is a semiconductor associated with the class, we take its  parameters
  if (semiconductor != NULL)
    {
      par = (dynamic_cast<ZbDDsemiconductor*>(semiconductor))->get_parameters();
    }

  if (min_name == "Gamma")
    {
      //--------------------------------------------------
      imass = Tensor2Sym(0);
      imass(1,1) = 1.0/par.m_G;
      imass(2,2) = 1.0/par.m_G;
      imass(3,3) = 1.0/par.m_G;

      edge = ( par.Ev + ((1.0/3.0) * par.delta ) + par.EgGamma)/Hartree;
      
    }
  //--------------------------------------------------------------------------------//
}

//=======================================================================//
SBZbCondBandBulkHamiltonian::SBZbCondBandBulkHamiltonian(ZbDDsemiconductor::ZbDDparameters& parameters):SBbulkHamiltonian()
{
  par = parameters;

  min_name = "Gamma";
  //---------------
  if (min_name == "Gamma")
    {
      imass = Tensor2Sym(0);
      imass(1,1) = 1.0/par.m_G;
      imass(2,2) = 1.0/par.m_G;
      imass(3,3) = 1.0/par.m_G;

      edge = ( par.Ev + ((1.0/3.0) * par.delta ) + par.EgGamma)/Hartree;

      kp_bands.resize(1,0);

    }
  
  
  //--------------
}


//===========================================================================//
void SBZbCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
 
  //now strain and potential
  if (min_name == "Gamma")
    {
       Hamiltonian[0][0].constant =  Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
	 + par.a_c * trace(strain_crystal )/Hartree;

      
    }
}
