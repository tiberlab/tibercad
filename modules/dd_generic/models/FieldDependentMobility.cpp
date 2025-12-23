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
 * \file FieldDependentMobility.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "FieldDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



void
FieldDependentMobility::read_database(void)
{
  const Database& db = get_database();
  db.set_section("mobility/field_dependent");

  std::vector<double> empty(2, 0);

  std::vector<double> data(empty);
  db.get("beta0", data, true);
  _beta =  get_carrier_type() == 'e' ? data[0] : data[1];

  data = empty;
  db.get("betaexp", data, true);
  _betaexp =  get_carrier_type() == 'e' ? data[0] : data[1];

  _vsat_formula = db.get("Vsat_Formula", _vsat_formula);

  if (_vsat_formula == 1)
  {
    data = empty;
    db.get("vsat0", data, true);
    _vsat0 =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("vsatexp", data, true);
    _vsat_b =  get_carrier_type() == 'e' ? data[0] : data[1];
  }
  else
  {
    data = empty;
    db.get("A_vsat", data, true);
    _vsat0 =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("B_vsat", data, true);
    _vsat_b =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("vsat_min", data, true);
    _vsat_min =  get_carrier_type() == 'e' ? data[0] : data[1];
  }

}


void
FieldDependentMobility::prepare_submodels(void)
{
  assert(_low_field_mob == NULL);

  ModelOptions opts;
  std::string low_field_model = get_option("low_field_model", "doping_dependent");
  if (get_options().has_submodel("low_field_model"))
  {
    opts = get_options().submodels_begin("low_field_model")->second;
    low_field_model = opts.get_name();
  }

  opts.set_option("type", low_field_model);


  opts.set_option("particle", get_option("particle", "electron"));

  create_submodel(_low_field_mob, "mobility", opts);
  _low_field_mob->set_carrier_type(get_carrier_type());
}


void
FieldDependentMobility::do_init(void)
{
  _beta = get_option("beta0", _beta);
  _betaexp = get_option("betaexp", _betaexp);
  _vsat_formula = get_option("Vsat_Formula", _vsat_formula);
  if (_vsat_formula == 1)
  {
    _vsat0 = get_option("vsat0", _vsat0);
    _vsat_b = get_option("vsatexp", _vsat_b);
  }
  else
  {
    _vsat0 = get_option("A_vsat", _vsat0);
    _vsat_b = get_option("B_vsat", _vsat_b);
    _vsat_min = get_option("vsat_min", _vsat_min);
  }

  _damping = get_option("damping_parameter", _damping);

  std::string force = get_option("driving_force", "grad_fermi");
  if (force == "efield")
    _force = EFIELD;
  else if (force == "grad_fermi")
    _force = GRADFERMI;
  else if (force == "field_parameter")
    _force = FIELDPARAM;
  else
  {
    std::string msg("FieldDependentMobility: Unknown driving force '");
    msg += force + "'.";
    throw InitFailedException(msg);
  }
}



double
FieldDependentMobility::get_mobility(void)
{
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;

  ID id = get_carrier();

  libMesh::RealGradient grad_fermi;
  get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);

  if (_force == GRADFERMI)
    E = grad_fermi.norm();
  else if (_force == EFIELD)
  {
    E = grad_fermi.norm();
    if (E > 1e-6)
      E = std::fabs(grad_fermi * get_driftdiffusionproperties().get_electric_field()) / E;
    else
      E = 0.0;
  }
  else if (_force == FIELDPARAM)
  {
    E = std::fabs(grad_fermi * get_driftdiffusionproperties().get_electric_field());
    E = std::sqrt(E);
  }

  double n = get_driftdiffusionproperties().get_q_density(id);

  double damp = n / (n + _damping);
  E *= damp;



  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);
  double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, beta);
  double mu = mu_lowfield * std::pow(tmp, -1.0/beta);

  return mu;
}




void
FieldDependentMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{

  ID id = get_carrier();

  dm.zero();
  double T = get_driftdiffusionproperties().get_lattice_temperature();

  libMesh::RealGradient grad_fermi;
  get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);
  double E = grad_fermi.norm();

  double n = get_driftdiffusionproperties().get_q_density(id);

  double damp = n / (n + _damping);

  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);

  if (0) // ((_force == FIELDPARAM))
  {
    libMesh::RealGradient grad_fermi;
    get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::sqrt(std::fabs(F));

    if (F > 1)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = (0.5 * dmu * sign * damp * damp / F) * grad_fermi;
    }
  }
  else if (0) //((_force == EFIELD) && (E > 1e-6))
  {
    libMesh::RealGradient grad_fermi;
    get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::fabs(F) / E;

    if (std::pow(mu_lowfield * F / vsat, beta) > 1e-3)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = dmu * sign * grad_fermi / E * damp;

    }

  }
}

void
FieldDependentMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{

  ID id = get_carrier();

  dm.zero();
  double T = get_driftdiffusionproperties().get_lattice_temperature();

  libMesh::RealGradient grad_fermi;
  get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);
  double E = grad_fermi.norm();

  double n = get_driftdiffusionproperties().get_q_density(id);

  double damp = n / (n + _damping);


  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);


  if ((_force == GRADFERMI) && (E > 1e-6))
  {
    E *= damp;

    double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, beta);
    double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
    dmu *= std::pow(mu_lowfield * E / vsat, beta) / (E * E);

    libMesh::RealGradient grad_fermi;
    get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);
    dm = dmu * grad_fermi * damp * damp;

  }
  else if (0) //((_force == FIELDPARAM))
  {
    libMesh::RealGradient grad_fermi;
    get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::sqrt(std::fabs(F));

    if (F > 1)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = get_driftdiffusionproperties().get_electric_field();
      dm *= 0.5 * dmu * sign * damp * damp / F;
    }
  }
  else if (0) //((_force == EFIELD) && (E > 1e-6))
  {
    libMesh::RealGradient grad_fermi;
    get_driftdiffusionproperties().get_grad_fermi(id, grad_fermi);

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::fabs(F) / E;

    //if (F > 1e-6)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = dmu * sign * get_driftdiffusionproperties().get_electric_field() / E * damp;

      dm -= dmu * F * grad_fermi / (E*E);
    }

  }
}
