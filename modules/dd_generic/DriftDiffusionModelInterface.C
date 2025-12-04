/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DriftDiffusionModelInterface.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "DriftDiffusionModelInterface.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/physics/Material.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/RuntimeException.h"

const double
DriftDiffusionModelInterface::T0 = Constants::k_B * 300.0;


DriftDiffusionProperties&
DriftDiffusionModelInterface::get_bulk_driftdiffusionproperties(void) const
{
  if (get_material() == NULL)
    throw RuntimeException("Drift-Diffusion model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<DriftDiffusionProperties&>(*
      get_material()->get_model(get_simulator_id()));

}

DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void) const
{
  if (get_owner() == NULL)
    throw RuntimeException("Drift-Diffusion model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<DriftDiffusionProperties&>(*
      get_owner()->get_model(get_simulator_id()));

}
