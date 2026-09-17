/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \file MTModel.C
 * \brief tiberCAD masstransport module implementation.
 *
 * \note This file is part of module masstransport.
 */


#include "MTModel.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/physics/Material.h"

using namespace std;

MTModel::~MTModel(void) = default;


TiberModelObject*
MTModel::_create(const ModelOptions& options)
{
  return new MTModel(options);
}


void
MTModel::_destroy(TiberModelObject* p)
{
  delete p;
}


MTModel*
MTModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  MTModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<MTModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    type = "bulk_" + type;
    pm = PhysicalModel::create<MTModel>(type, mat, options);
  }

  return(pm);
}



void
MTModel::do_init(void)
{

  // we read it in g/m^3/Pa
  _solubility = get_option("solubility", _solubility);

  // we read it in m^2/S
  double d0 = get_option("diffusivity_0", 3.21e-4);
  double Ea = get_option("diff_act_en",0.415); //eV
  const double kb = 8.617e-5; //eV/K
  _cell_temp = get_option("cell_temperature", SimulationOptions::temperature);

  //Arrhenius D = D0 exp(-Ea/kbT)
  _diffusivity = d0*exp(-Ea/(kb*_cell_temp));

}


void
MTModel::calculate(const Elem* elem, const Point& point)
{
 
}

void
MTModel::prepare_submodels(void)
{
}



