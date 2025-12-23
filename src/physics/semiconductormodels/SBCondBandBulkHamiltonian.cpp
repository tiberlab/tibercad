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
 * \file SBCondBandBulkHamiltonian.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/semiconductormodels/SBCondBandBulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/Semiconductor.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Messages.h"

using namespace std;

SBCondBandBulkHamiltonian::~SBCondBandBulkHamiltonian(void)
{
}



void
SBCondBandBulkHamiltonian::prepare_submodels(void)
{
  assert(semiconductor == NULL);

  ModelOptions opt =  get_options();
  opt.delete_all_submodels();
  semiconductor = Semiconductor::create(get_material(), opt);
  add_submodel("semiconductor", semiconductor);
}



void
SBCondBandBulkHamiltonian::do_init(void)
{
  SBbulkHamiltonian::do_init();

  kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(0,0));
  // NOTE: this duplicates the band, which makes OpticsKP use two 
  // conduction bands, giving the same result as an 8x8 kp. The other
  // parts of the code are not affected by this doubling, as they do not refer
  // to kp_bands_map
  kp_bands_map.insert(std::make_pair(1,0));

  calculate_for_init();

  get_option("model", "");
}




