using namespace std;
#include "SBWzCondBandBulkHamiltonian.h"
#include "getpot.h"
#include "Alloy.h"

//=======================================================================//
void SBWzCondBandBulkHamiltonian::do_init()
{

  SBCondBandBulkHamiltonian::do_init();

  WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor)) -> get_parameters();

  wz_par = &par;


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

  calculate_Hamiltonian_gen();
  
}





//===========================================================================//
void SBWzCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
 
  //now strain and potential

  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
    + ((strain_crystal(1,1) + strain_crystal(2,2))* (wz_par->a_x) + (wz_par->a_z) *  strain_crystal(3,3))/Hartree;
    
}

//============================================================================//



