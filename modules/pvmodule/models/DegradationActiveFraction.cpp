/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file DegradationActiveFraction.C
 * \brief tiberCAD pvmodule module implementation.
 *
 * \note This file is part of module pvmodule.
 */

#include "DegradationActiveFraction.h"
#include "tibercad/module/SimulationInterface.h"

#include "tibercad/module/TiberModule.h"


DegradationActiveFraction::DegradationActiveFraction(const ModelOptions& options)
  : DegradationModel(options)
{
}

DegradationActiveFraction*
DegradationActiveFraction::create(const ModelOptions& options)
{
  return new DegradationActiveFraction(options);
}

void
DegradationActiveFraction::do_init(void)
{

  if (get_options().has_submodel("photocurrent"))
  {
    ModelOptions& opts = get_options().submodels_begin("photocurrent")->second;
    _exponent_ph = opts.get_option("exponent", _exponent_ph);
  }

  if (get_options().has_submodel("series_resistance"))
  {
    ModelOptions& opts = get_options().submodels_begin("series_resistance")->second;
    _exponent_rs = opts.get_option("exponent", _exponent_rs);
  }
 
  std::string degradation_kinetic_sim = get_option("psk_fraction", "");
  _kinetic_model = SimulationInterface::find_solution_provider(degradation_kinetic_sim, "PerovskiteFraction");
}

void
DegradationActiveFraction::do_degrade_params(const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  DegradationModel::Parameters& params) const
{
  double psk_f = 1;

  if (_kinetic_model.first != nullptr)
  { 
    if (_kinetic_model.first->get_solution(elem, _kinetic_model.second, psk_f, p, false))
    {
     // factor for photocurrent
     double ph_fac = std::pow(psk_f , _exponent_ph);
     params.double_params[2] *= ph_fac; //J = J0 x f^n

     // factor for series resistance
     double rs_fac = std::pow(psk_f , _exponent_rs);
     params.double_params[0] /= rs_fac; //Rs = Rs,0 / f^n  
    }
  }
}




