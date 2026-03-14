/*  
 * This file is part of the tiberCAD module dd_generic.
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
 * \file FieldAssistedMobility.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "FieldAssistedMobility.h"
#include "DriftDiffusionProperties.h"

#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



void
FieldAssistedMobility::read_database(void)
{
  const Database& db = get_database();
  db.set_section("mobility/field_assisted");

  std::vector<double> mu0(2, _mu0);

  db.get("mu0", mu0);
  _mu0 =  get_carrier_type() == 'e' ? mu0[0] : mu0[1];

  std::vector<double> E0(2, _E0);
  db.get("E0", E0);
  _E0 =  get_carrier_type() == 'e' ? E0[0] : E0[1];
}



void
FieldAssistedMobility::do_init(void)
{
  get_parameter("mu0", _mu0);
  get_parameter("E0", _E0);

//  std::string force = get_option("driving_force", "efield");
//  if (force == "efield")
//    _force = EFIELD;
//  else if (force == "grad_fermi")
//    _force = GRADFERMI;
//  else
//  {
//    std::string msg("FieldAssistedMobility: Unknown driving force '");
//    msg += force + "'.";
//    throw InitFailedException(msg);
//  }

}



double
FieldAssistedMobility::get_mobility(void)
{
  //double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;

  /*
  const libMesh::RealGradient& grad_fermi = (get_carrier_type() == 'e') ?
    get_driftdiffusionproperties().get_grad_fermi_e() :
    get_driftdiffusionproperties().get_grad_fermi_h();

  E = grad_fermi.size();

  if ((_force == EFIELD) && (E > 1e-6))
    E = grad_fermi * get_driftdiffusionproperties().get_electric_field() / E;
  */

  E = get_driftdiffusionproperties().get_electric_field().norm();

  double arg = std::sqrt(E / _E0);

  double mu = _mu0 * std::exp(arg);

  return mu;
}



void
FieldAssistedMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const libMesh::RealGradient& E = dd.get_electric_field();
  double mu  = get_mobility();

  double dmu = (E.norm() > 0) ? mu / (2.0 * E.norm() * sqrt(E.norm() * _E0)) : 0.0;
  dm = -E * dmu;
}

void
FieldAssistedMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  //std::cout<<"#############################____##############################"<<std::endl;
  dm.zero();
  //double T = get_driftdiffusionproperties().get_lattice_temperature();
  //double E = 0.0;
  //if (get_carrier_type() == 'e')
  //  E = get_driftdiffusionproperties().get_grad_fermi_e().size();
  //else
  //  E = get_driftdiffusionproperties().get_grad_fermi_h().size();
/*
  if ((_force == EFIELD) && (E > 1.0))
  {

    double mu = get_mobility();
    E = get_driftdiffusionproperties().get_electric_field().size();
    double tmp = 0.5 / (std::sqrt(E * _E0) * E);


    dmu *= mu * tmp * get_driftdiffusionproperties().get_electric_field();

    if (get_carrier_type() == 'e')
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_e();
    else
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_h();

  }
*/
}

