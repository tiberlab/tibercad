/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file AugerRecombination.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "AugerRecombination.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/io/Database.h"

#include "tibercad/module/TiberModule.h"


void
AugerRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/auger");

  std::vector<double> data(2, 0);
  db.get("A", data, true);
  _An = data[0];
  _Ap = data[1];

  data = std::vector<double>(2, 0);
  db.get("B", data);
  _Bn = data[0];
  _Bp = data[1];

  data = std::vector<double>(2, 0);
  db.get("C", data);
  _Cn = data[0];
  _Cp = data[1];

  data = std::vector<double>(2, 0);
  db.get("H", data);
  _Hn = data[0];
  _Hp = data[1];

  data = std::vector<double>(2, 0);
  db.get("N0", data);
  _N0n = data[0];
  _N0p = data[1];

}



void
AugerRecombination::do_init(void)
{

  if (has_parameter("Cn"))
  {
    get_parameter("Cn", _Cn);
    _fixed_Cn = true;
  }
  if (has_parameter("Cp"))
  {
    get_parameter("Cp", _Cp);
    _fixed_Cp = true;
  }
}


inline
double
AugerRecombination::get_Cn(void)
{
  double Cn = _Cn;
  if (!_fixed_Cn)
  {
    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
    double T = dd.get_lattice_temperature() / T0;
    double tmp = _An + _Bn * T + _Cn * T * T;
    Cn = tmp * std::exp(-dd.get_electron_density() / _N0n);
  }
  return Cn;
}

inline
double
AugerRecombination::get_Cp(void)
{
  double Cp = _Cp;
  if (!_fixed_Cp)
  {
    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
    double T = dd.get_lattice_temperature() / T0;
    double tmp = _Ap + _Bp * T + _Cp * T * T;
    Cp = tmp * std::exp(-dd.get_electron_density() / _N0p);
  }
  return Cp;
}



void
AugerRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature();
  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double np = n * p;

  long double g = 1.0 - exp((Efp - Efn) / T);

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * n;
  long double B = Cp * p;
  recomb_e = recomb_h = (A + B) * np * g;
}



void
AugerRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature();
  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double np = n * p;

  long double e = exp((Efp - Efn) / T);
  long double g = 1 - e;

  double Cn = get_Cn();
  double Cp = get_Cp();

  long double A = Cn * n;
  long double B = Cp * p;

  long double dRedn = (2*A + B) * p * g;
  long double dRedp = (2*B + A) * n * g;
  long double dRedEfn = -(A + B) * np * e / T;
  long double dRedEfp = -dRedEfn;

  recomb_e[0] = recomb_h[0] = dRedn;
  recomb_e[1] = recomb_h[1] = dRedp;
  recomb_e[2] = recomb_h[2] = dRedEfn;
  recomb_e[3] = recomb_h[3] = dRedEfp;
}




