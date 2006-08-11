using namespace std;
#include "SBWzCondBandBulkHamiltonian.h"
#include "getpot.h"
#include "Alloy.h"




//===========================================================//
SBWzCondBandBulkHamiltonian::SBWzCondBandBulkHamiltonian( ):SBbulkHamiltonian()
{
  kp_bands.resize(1,0);
  kp_bands_map.insert(make_pair(0,0));
}

SBWzCondBandBulkHamiltonian::SBWzCondBandBulkHamiltonian(const SBWzCondBandBulkHamiltonian &  model):SBbulkHamiltonian(model)
{ 
  kp_bands.resize(1,0);
  par = model.par;
  calculate_edge_and_mass();
  kp_bands_map.insert(make_pair(0,0));
}

//=======================================================================//
void  SBWzCondBandBulkHamiltonian::calculate_edge_and_mass()
{
  //--------------------------------------------------
  //if there is a semiconductor associated with the class, we take its  parameters
  if (semiconductor != NULL)
    {
      par = (dynamic_cast<WzDDsemiconductor*>(semiconductor))->get_parameters();
    }


  //--------------------------------------------------
  imass = Tensor2Sym(0);
  imass(1,1) = 1.0/par.m_c_xx;
  imass(2,2) = 1.0/par.m_c_xx;
  imass(3,3) = 1.0/par.m_c_zz;

  //--------------------------------
  double d1 =  par.delta_cr;
  double d2 =  par.delta_s;
  double d3 =  d2;
  
  double E1 = d1 + d2;

  double E2 = (d1 - d2)/2.0 + sqrt( (d1-  d2/2.0)*( d1- d2/2.0) + 2.0 * d3 * d3 );

  double Ev_top;

  if (E1 > E2)
    Ev_top = par.Ev + E1;
  else
    Ev_top = par.Ev + E2;
  
  //--------------------------------------------------------------------------------/

  edge = (Ev_top + par.EgGamma)/Hartree;

  //--------------------------------------------------------------------------------//
}

//=======================================================================//
SBWzCondBandBulkHamiltonian::SBWzCondBandBulkHamiltonian(WzDDsemiconductor::WzDDparameters& parameters):SBbulkHamiltonian()
{
  kp_bands.resize(1,0);
  kp_bands_map.insert(make_pair(0,0));
  par = parameters;
  calculate_edge_and_mass();
  
  
}


//===========================================================================//
void SBWzCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
 
  //now strain and potential

  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
    + ((strain_crystal(1,1) + strain_crystal(2,2))* par.a_x + par.a_z *  strain_crystal(3,3))/Hartree;
    
}

//============================================================================//



