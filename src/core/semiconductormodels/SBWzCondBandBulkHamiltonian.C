// $Id$

#include "SBWzCondBandBulkHamiltonian.h"
#include "Material.h"
#include "Constants.h"

using namespace std;
using namespace Constants;



//===========================================================================//
void SBWzCondBandBulkHamiltonian::calculate_for_init( )
{
  const WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor)) -> get_parameters();

  wz_par = &par;


   //--------------------------------------------------
  imass = Tensor2Sym(0);
  imass(1,1) = 1.0/par.m_c_xx;
  imass(2,2) = 1.0/par.m_c_xx;
  imass(3,3) = 1.0/par.m_c_zz;


  //--------------------------------------------------------------------------------/

  edge = (par.Ev + par.EgGamma) / Hartree;

  //--------------------------------------------------------------------------------//

  calculate_Hamiltonian_gen();

  calculate_Hamiltonian_k_par();

}


//===========================================================================//
void SBWzCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{

  //now strain and potential

  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
    + ((strain_crystal(1,1) + strain_crystal(2,2))* (wz_par->a_x) + (wz_par->a_z) *  strain_crystal(3,3))/Hartree;

}

//============================================================================//



