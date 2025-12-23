/*  
 * This file is part of the tiberCAD module dssc.
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
 * \file DSSCContact.C
 * \brief tiberCAD dssc module implementation.
 *
 * \note This file is part of module dssc.
 */


#include "DSSCContact.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/physics/Constants.h"

using namespace std;

bool
DSSCContact::_open_circuit = true;




DSSCContact::DSSCContact(const ModelOptions& options)
  // open circuit value
  : PhysicalModel(options),
    //_boundary_value(1e10),
    _bias(0.0),
    _res(1e10),
    _cathode(false),
    _gate(false),
    _j0(0.1),
    _beta(1.0),
    _current(0.0),
    _barrier(0.0),
    _kinetic_rate(1e4),
    _kinetic_rate_I3(0),
    _kinetic_rate_I(0)
{
}



DSSCContact*
DSSCContact::create(const MaterialBoundary* boundary,
    const ModelOptions& options)
{
  DSSCContact* ct = NULL;

  string name = options.get_option("type", "");

  ct = new DSSCContact(options);

  if (ct != NULL)
  {
    ct->set_options(options);

    if (name == "Pt")
      ct->is_cathode() = true;
    else if (name == "Gate")
      ct->is_gate() = true;
    else if (name != "ohmic")
      throw ModelErrorException(
          "DSSC: No such boundary model: " + name);
  }
  else
    throw ModelErrorException("Could not create DSSC boundary model");

  return ct;
}



void
DSSCContact::do_init(void)
{
  get_parameter("load", _res);
  get_parameter("bias", _bias);
  get_parameter("barrier", _barrier);
  get_parameter("kinetic_rate", _kinetic_rate);
  get_parameter("kinetic_rate_I", _kinetic_rate_I);
  get_parameter("kinetic_rate_I3", _kinetic_rate_I3);
  get_parameter("Ex_curr", _j0);
  //get_parameter("beta", _beta);
}



void
DSSCContact::calculate_current(double I, double I3)
{
  //double kT = Constants::k_B * SimulationOptions::T;
  //double upt = 0.0;
  //double A = sqrt(I3 * _Ioc / (I * _I3oc)) * exp((1 - _beta) * upt / kT);
  //double B = I / _Ioc * exp(-_beta * upt / kT);

  //_current = _j0 * (A - B);
  //_current = get_potential();
}
