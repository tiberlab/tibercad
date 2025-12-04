/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file SBWzCondBandBulkHamiltonian.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/semiconductormodels/SBWzCondBandBulkHamiltonian.h"
#include "tibercad/physics/Material.h"
#include "tibercad/physics/Constants.h"

using namespace std;
using namespace Constants;



//===========================================================================//
void SBWzCondBandBulkHamiltonian::calculate_for_init( )
{
  const WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor)) -> get_parameters();

  wz_par = &par;


   //--------------------------------------------------
  imass = Tensor2(0);
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
void SBWzCondBandBulkHamiltonian::apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential)
{

  //now strain and potential

  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
    + ((strain_crystal(1,1) + strain_crystal(2,2))* (wz_par->a_x) + (wz_par->a_z) *  strain_crystal(3,3))/Hartree;

}

//============================================================================//



