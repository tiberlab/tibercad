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
 * \file NonlinearPolarization.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "NonlinearPolarization.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"
#include "tibercad/math/TensorOperators.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


void
NonlinearPolarization::do_init(void)
{
  PolarizationModel::do_init();

  std::string sim_name = "";
  get_parameter("strain_simulation", sim_name);
  _strain.set_simulation(sim_name);

}


void
NonlinearPolarization::read_database(void)
{
  // get piezoelectric_coefficients
  const Database& db = get_database();
      
  if (get_material()->get_structure() == "wz")
  {
    db.set_section("polarization/PRB88");
    _e33 = db.get("e33", 0.0);
    _e31 = db.get("e31", 0.0);
    _e15 = db.get("e15", 0.0);

    _2nd_order_coeff.resize(8, 0.0);
    db.get("2nd_order_coefficients", _2nd_order_coeff);

    _2nd_order_coeff[0] = db.get("B115", _2nd_order_coeff[0]);
    _2nd_order_coeff[1] = db.get("B125", _2nd_order_coeff[1]);
    _2nd_order_coeff[2] = db.get("B135", _2nd_order_coeff[2]);
    _2nd_order_coeff[3] = db.get("B311", _2nd_order_coeff[3]);
    _2nd_order_coeff[4] = db.get("B312", _2nd_order_coeff[4]);
    _2nd_order_coeff[5] = db.get("B313", _2nd_order_coeff[5]);
    _2nd_order_coeff[6] = db.get("B333", _2nd_order_coeff[6]);
    _2nd_order_coeff[7] = db.get("B344", _2nd_order_coeff[7]);

    _Psp = db.get("Pz", _Psp);
  } 
  else if (get_material()->get_structure() == "zb")
  {
    db.set_section("piezoelectricity");
    _e14 = db.get("e14", 0.0);

    _2nd_order_coeff.resize(3, 0.0);
    db.get("2nd_order_coefficients", _2nd_order_coeff);

    _2nd_order_coeff[0] = db.get("B114", _2nd_order_coeff[0]);
    _2nd_order_coeff[1] = db.get("B124", _2nd_order_coeff[1]);
    _2nd_order_coeff[2] = db.get("B156", _2nd_order_coeff[2]);
  }

}


void
NonlinearPolarization::do_calculate(const Elem* elem, const Point& point)
{  

  RealVectorValue polarization(0);

  Tensor2& strain = get_strain();
  _strain.get_crystal_strain(elem, point, strain);

  // strain in Voigt notation
  double e1 = strain(1,1);
  double e2 = strain(2,2);
  double e3 = strain(3,3);
  double e4 = 2*strain(3,2);
  double e5 = 2*strain(3,1);
  double e6 = 2*strain(2,1);

  // compute polarization
  if (get_material()->get_structure() == "wz")
  {

    // a = 0, b = 1, .. h = 7
    double B125 = _2nd_order_coeff[1];
    double B115 = _2nd_order_coeff[0];
    double B135 = _2nd_order_coeff[2];
    double B146 = 0.5 * (_2nd_order_coeff[0] - _2nd_order_coeff[1]);
    double B214 = B125;
    double B224 = B115;
    double B234 = B135;
    double B256 = B146;
    double B311 = _2nd_order_coeff[3];
    double B322 = B311;
    double B333 = _2nd_order_coeff[6];
    double B344 = _2nd_order_coeff[7];
    double B355 = B344;
    double B366 = 0.5 * (_2nd_order_coeff[3] - _2nd_order_coeff[4]);
    double B312 = _2nd_order_coeff[4];
    double B313 = _2nd_order_coeff[5];
    double B323 = B313;

    polarization(0) = _e15 * e5;
    polarization(1) = _e15 * e4;
    polarization(2) = _e31 * e1 + _e31 * e2 + _e33 * e3;

    polarization(0) += 0.5 * (B125 * e2 * e5 + B115 * e1 * e5 +
                              B135 * e3 * e5 + B146 * e4 * e6);
    polarization(1) += 0.5 * (B214 * e1 * e4 + B224 * e2 * e4 +
                              B234 * e3 * e4 + B256 * e5 * e6);
    polarization(2) += 0.5 * (B311 * e1 * e1 + B322 * e2 * e2 + B333 * e3 * e3 +
                              B344 * e4 * e4 + B355 * e5 * e5 + B366 * e6 * e6 +
                              B312 * e1 * e2 + B313 * e1 * e3 + B323 * e2 * e3);


    polarization(2) += _Psp;
  }
  else
  {
    double B114 = _2nd_order_coeff[0];
    double B124 = _2nd_order_coeff[1];
    double B156 = _2nd_order_coeff[2];

    polarization(0) = _e14 * e4;
    polarization(1) = _e14 * e5;
    polarization(2) = _e14 * e6;

    polarization(0) += B114 * e1 * e4 + B124 * e4 * (e2 + e3) + B156 * e5 * e6;
    polarization(1) += B114 * e2 * e5 + B124 * e5 * (e3 + e1) + B156 * e4 * e6;
    polarization(2) += B114 * e3 * e6 + B124 * e6 * (e1 + e2) + B156 * e4 * e5;
  }

  set_polarization(polarization);

  // rotate polarization to calculation system
  rotate();

}
