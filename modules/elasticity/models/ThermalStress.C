/*  
 * This file is part of the tiberCAD module elasticity.
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
 * \file ThermalStress.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "ThermalStress.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/math/TensorOperators.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


ThermalStress::ThermalStress(const ModelOptions& options) :
    BodyForceModel(options)
{
}

ThermalStress::~ThermalStress(void)
{

}


void
ThermalStress::do_init(void)
{

  libMesh::RealGradient body_force(0);
  
 //Get reference lattice
  _ref_temp = get_option("reference_temperature", SimulationOptions::temperature);
 
  get_parameter("thermal_expansion_coefficient", _alpha);

 
  std::string temp_sim = get_option("thermal_simulation", "");

  _temp.set_simulation(temp_sim);
}



void
ThermalStress::read_database(void)
{

  const Database& db = get_database();
  db.set_section("lattice");

  vector<double> alpha;
  db.get("thermal_coefficient", alpha);

  switch (alpha.size())
  {
    case 1:
      _alpha(0) = _alpha(1) = _alpha(2) = alpha[0];
      break;
    case 2:
      _alpha(0) = _alpha(1) = alpha[0];
      _alpha(2) = alpha[1];
      break;
    case 3:
      _alpha(0) = alpha[0];
      _alpha(1) = alpha[1];
      _alpha(2) = alpha[2];
      break;
    default:
      break;
  }

}

 

void
ThermalStress::calculate(const libMesh::Elem* elem, const libMesh::Point& point)
{
  // get temperature
  double deltaT = _temp.get_temperature(elem, point, true) - _ref_temp;

  const libMesh::RealTensor& rotm = get_material()->get_rotation_matrix();


  // compute thermally induced strain
  libMesh::RealTensor strain(0);
  strain(0,0) = -_alpha(0) * deltaT;
  strain(1,1) = -_alpha(1) * deltaT;
  strain(2,2) = -_alpha(2) * deltaT;

  // rotate
  strain = rotm * (strain * rotm.transpose());

  set_strain_source(strain);
}

