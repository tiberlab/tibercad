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
 * \file SBuserHamiltonian.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/semiconductormodels/SBuserHamiltonian.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/math/TensorOperators.h"


using namespace std;

//======================================================================//
SBuserHamiltonian::~SBuserHamiltonian()
{
  
}

//=======================================================================//
SBuserHamiltonian::SBuserHamiltonian(const ModelOptions& options)
  : SBbulkHamiltonian(options)
{
  edge = 0.0;
 
  
  imass = Tensor2(1);
}



//======================================================================//


void SBuserHamiltonian::do_init( )
{


  EFAbulkHamiltonian::do_init();

  // for now this is automatically a valence band
  if (get_option("particle", "") == "hl")
  {
    // band edge
    const Database& db = get_database();
    db.set_section("valenceband");
    edge = db.get("E_v", 0.0);

    // mass
    const string& band = get_option("band", "hh");
  }

  // override of band edge
  edge = get_option("band_edge", edge);
  edge /= Constants::Hartree;

  // one degenerate band
  kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(2,0));
  // NOTE: this duplicates the band, which makes OpticsKP use two
  // identical bands with opposite spin
  // BUT: maybe optics cannot do anything with this single band anyways?
  kp_bands_map.insert(std::make_pair(3,0));

  double mass = get_option("effective_mass", 1.0);

  if (mass == 0.0) throw InitFailedException("User-defined Hamiltonian: zero mass");

  imass = -(1.0 / mass) * Tensor2(1.0);
  

  calculate_Hamiltonian_gen();
  
}

//======================================================================//
