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
 * \file ConstantDOS.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */

#include "ConstantDOS.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/physics/Material.h"
#include "tibercad/math/TiberMath.h"
#include "tibercad/io/Database.h"
#include "tibercad/base/InitFailedException.h"

#include "tibercad/module/TiberModule.h"

using namespace std;

ConstantDOS::ConstantDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _Ewidth(0.0), _N0(1.0)
{

}

void
ConstantDOS::read_database(void)
{
  // when reading from the database, we use the same data
  // as for kp
  const Database& db = get_database();

  if (get_particle() == 'e')
  {
    // TODO: should bowing be applied to Eg(T) or Eg(0) ?

    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);

    db.set_section("bandgap");
    double bandgap = db.get("Eg_G", 1e3);

    reference_energy()[0] += bandgap;
  }
  else if (get_particle() == 'h')
  {
    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);
  }
}

void
ConstantDOS::do_init(void)
{
  if (has_parameter("level"))
    get_parameter("level", reference_energy());

  get_parameter("Ewidth", _Ewidth);
  get_parameter("N0", _N0);

  effective_mass()[0] = 1.0;

  effective_dos() = _N0;
  total_state_density() = _N0;
}

void
ConstantDOS::calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
    double kT, double , const Elem* , const Point& ) const
{

  double dens, der, der2, Emin, Emax, exp_num, exp_den, ref_en;
  //cout<<"ref_energy "<<reference_energy()<<endl;

  ref_en = get_reference_energy()[0];
  if (get_particle() == 'h') ref_en *= -1.0;

  if (get_particle() == 'e')
  {
    Emin = ref_en;
    Emax = Emin + _Ewidth; 
  }
  else
  {
    Emax = ref_en;
    Emin = Emax + _Ewidth;
  }
  
  double zmin = (Ef - Emin - Epot) / kT;
  double zmax = (Ef - Emax - Epot) / kT;
  
  if (get_particle() == 'e')
  {
    exp_num = exp(zmin);
    exp_den = exp(zmax);
  }
  else
  {
    exp_num = exp(zmax);
    exp_den = exp(zmin);
  }

  double C = _N0 / _Ewidth;

  dens = kT * C * log((exp_num + 1)/(exp_den +1));
  if (result.size() > 1)
    der = C * ( exp_num / (exp_num + 1) - exp_den / (exp_den + 1));
  if (result.size() > 2)
    der2 = C/kT * (std::pow(exp_num/(exp_num + 1),2) - std::pow(exp_den/(exp_den + 1),2));

  
  result[0] = dens;
  if (result.size() > 1)
    result[1] = der;
  if (result.size() > 2)
    result[2] = der2;

}
