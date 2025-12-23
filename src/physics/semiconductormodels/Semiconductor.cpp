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
 * \file Semiconductor.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/semiconductormodels/Semiconductor.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Alloy.h"
#include "tibercad/base/SimulationOptions.h"



using namespace std;

Semiconductor::Semiconductor(const ModelOptions& options)
 : PhysicalModel(options)
{
  modelA = NULL;

  modelB = NULL;

  //database is read before do_init() and this parameter must be ready.
  _kp_model = get_option("model", "6x6");
  _consider_temperature = get_option("consider_temperature",  true );
  _consider_temperature = get_option("temperature_scaling", _consider_temperature );
  _couple_bands = get_option("cb_vb_coupling", true );
}


Semiconductor* Semiconductor::create(const Material* mat,  const ModelOptions& options)
{
  std::string structure = mat->get_structure();
  return PhysicalModel::create<Semiconductor>("semicond_" + structure, mat, options);
}


//--------------------------------------------------------------------------------------------//
void Semiconductor::do_init ()
{

  int verbose = SimulationOptions::verbose ();

  get_option("particle","");
  get_option("particle","");
  get_option("kpVVtermSymmetric","");
  get_option("kpCVtermSymmetric","");

  // the temperature simulation
  string temp_simul = get_option("thermal_simulation", "");

  temp_interface.set_simulation(temp_simul);


  _spurious = get_option("spurious","Chuang");

  _temperature = 0.0;
  if (_consider_temperature)
  { 
      _temperature = get_option("temperature", SimulationOptions::T);
  }
  

  /*

  if  (verbose > 0)
    if (_consider_temperature)
      std::cout << "Semiconductor: band gap depends on temperature\n";
    else
      std::cout << "Semiconductor: band gap does not depend on temperature\n";

  if  (verbose > 0)
    if (temp_simul != "")
      std::cout << "Semiconductor: temparature is taken from the simulation " << temp_simul << "\n";


  */
}





//---------------------------------------------------------------------------------------------//
inline
void Semiconductor::do_init_alloy (const PhysicalModel *comp_A, const PhysicalModel *comp_B, double xa)
{

  modelA = dynamic_cast<const Semiconductor* >(comp_A);

  modelB = dynamic_cast<const Semiconductor* >(comp_B);

  _xa = xa;

}
