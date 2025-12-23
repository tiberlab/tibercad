/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file ExcitonGeneration.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "ExcitonGeneration.h"

#include "DriftDiffusionProperties.h"

#include "tibercad/module/TiberModule.h"



void
ExcitonGeneration::do_init(void)
{
  C_ = get_option("C", 1e-10);
}

void
ExcitonGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e = recomb_h = C_ * n * p;
}

void
ExcitonGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = C_ * p;
  recomb_e[1] = recomb_h[1] = C_ * n;
}
