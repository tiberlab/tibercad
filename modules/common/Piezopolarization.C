/*  
 * This file is part of the tiberCAD module common.
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
 * \file Piezopolarization.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "Piezopolarization.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"
#include "tibercad/math/TensorOperators.h"


#include "tibercad/module/TiberModule.h"


using namespace std;


void
Piezopolarization::do_init(void)
{
  PolarizationModel::do_init();

  std::string sim_name = "";
  get_parameter("strain_simulation", sim_name);
  _strain = SimulationInterface::find_solution_provider(sim_name, "Strain");

}


void
Piezopolarization::read_database(void)
{
  // get piezoelectric_coefficients
  const Database& db = get_database();
  db.set_section("piezoelectricity");
      
  if (get_material()->get_structure() == "wz")
  {
    _e33 = db.get("e33", 0.0);
    _e31 = db.get("e31", 0.0);
    _e15 = db.get("e15", 0.0);
  } 
  else if (get_material()->get_structure() == "zb")
    _e14 = db.get("e14", 0.0);

}


void
Piezopolarization::do_calculate(const Elem* elem, const Point& point)
{  

  RealVectorValue polarization(0);

  libMesh::RealTensor strain;
  get_strain(strain);

  if (_strain.is_valid())
  {
    std::vector<Point> p(1);
    p[0] = point;
    std::vector<double> values(6);

    if (_strain.simulation()->get_solution(elem, _strain.id(), values, p))
    {
      strain(0,0) = values[0];
      strain(1,1) = values[1];
      strain(2,2) = values[2];
      strain(1,0) = strain(0,1) = values[3];
      strain(2,1) = strain(1,2) = values[4];
      strain(2,0) = strain(0,2) = values[5];
    }

    const Material* mat = get_material();
    const libMesh::RealTensor& rotm = mat->get_rotation_matrix();
    strain = rotm.transpose() * strain * rotm;
  }
  

  // compute polarization
  if (get_material()->get_structure() == "wz")
  {
    polarization(0) = 2.0 * _e15 * strain(2,0);
    polarization(1) = 2.0 * _e15 * strain(2,1);
    polarization(2) = _e31 * strain(0,0) + _e31 * strain(1,1) + _e33 * strain(2,2);
  }
  else
  {
    polarization(0) = 2.0 * _e14 * strain(2,1);
    polarization(1) = 2.0 * _e14 * strain(2,0);
    polarization(2) = 2.0 * _e14 * strain(1,0);
  }
  set_polarization(polarization);

  // rotate polarization to calculation system
  rotate();

}
