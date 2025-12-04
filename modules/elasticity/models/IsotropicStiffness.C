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
 * \file IsotropicStiffness.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "IsotropicStiffness.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


IsotropicStiffness::IsotropicStiffness(const ModelOptions& options):StiffnessModel(options)
{
}


// void IsotropicStiffness::read_database( )
// {

//   Database& db = get_database();
//   db.set_section("stiffness/isotropic");

//   _young = db.get("young", 0.0, false);
//   _poisson = db.get("poisson", 0.0, false);
 

// }


void
IsotropicStiffness::do_init(void)
{

  _young = 0.0;
  _poisson = 0.0;

  get_parameter("young",_young);
  get_parameter("poisson",_poisson);


  Tensor4DSym stiffness(0);

  double A = _young / (1 + _poisson) / (1 - 2 * _poisson);


  stiffness(1,1,1,1) = 1 - _poisson;
  stiffness(2,2,2,2) = 1 - _poisson;
  stiffness(3,3,3,3) = 1 - _poisson;
  stiffness(2,2,1,1) = _poisson;
  stiffness(3,3,1,1) = _poisson;
  stiffness(3,3,2,2) = _poisson;

  stiffness(3,2,3,2) = (1 - 2.0 * _poisson)/2.0;
  stiffness(3,1,3,1) = (1 - 2.0 * _poisson)/2.0;
  stiffness(2,1,2,1) = (1 - 2.0 * _poisson)/2.0;

  stiffness *= A;

  set_stiffness_constant(stiffness);

  //rotate();
  
}




 
