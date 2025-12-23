/*  
 * This file is part of the tiberCAD module thermal.
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
 * \file ThermalModel.C
 * \brief tiberCAD thermal module implementation.
 *
 * \note This file is part of module thermal.
 */


#include "ThermalModel.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "tibercad/io/Messages.h"
#include "tibercad/physics/misc/HeatSourceModel.h"
#include "tibercad/physics/misc/ThermalConductivityModel.h"


using namespace std;



ThermalModel::ThermalModel(const ModelOptions& options)
  : PhysicalModel(options),
    _tcm(NULL),
    _kappa(0),
    _heat_source(0)
{

}


ThermalModel*
ThermalModel::create(const Material* mat, const ModelOptions& options)
{
  return PhysicalModel::create<ThermalModel>(_create,_destroy, mat, options);
}

void
ThermalModel::prepare_submodels(void)
{
  
  ModelOptions opts;
  opts.set_option("type","constant");
  create_submodel(_tcm, "thermal_conductivity", opts);
  create_submodels(_hsm, "heat_source");

}

void
ThermalModel::do_init(void)
{
  _kappa = _tcm->get_thermal_conductivity();
}



void
ThermalModel::calculate(const Elem* elem, const Point& point, double temperature)
{
  //Heat Source
  _heat_source = 0.0;
 
  for (ID n = 0 ; n <_hsm.size() ; n++)
  {
    _hsm[n]->calculate(elem,point);
    _heat_source +=  _hsm[n]->get_heat_source();
  }

  _tcm->calculate(elem, point, temperature);
  _kappa = _tcm->get_thermal_conductivity();

}


    //! Print some useful information
void 
ThermalModel::do_print_info(void)
{
  Messages::info("Thermal conductivity:");
  ostringstream os;
  os <<"  Kxx: "<<_kappa(0,0)<<" W/(m K)\n";
  os <<"  Kyy: "<<_kappa(1,1)<<" W/(m K)\n";
  os <<"  Kzz: "<<_kappa(2,2)<<" W/(m K)";
  Messages::info(os.str());
  Messages::newline();
}
