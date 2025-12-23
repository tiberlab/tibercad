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
 * \file BandProperties.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "BandProperties.h"
#include "tibercad/physics/misc/ParticleDensity.h"
#include "tibercad/physics/misc/DensityOfStates.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/base/ModelOptions.h"
#include "tibercad/io/Messages.h"
#include "tibercad/io/Database.h"

#include "tibercad/module/TiberModule.h"

#include <sstream>
#include <algorithm>

using namespace std;



BandProperties::BandProperties(const ModelOptions& options) :
    DriftDiffusionModelInterface(options),
    _particle('e'),
    _dos_factor(pow(2.0 * M_PI *
        Constants::me / (Constants::h * Constants::h) *
        Constants::e, 1.5) / 1e6)
{
  string particle(get_option("particle", "-"));

  if (particle == string("el") || particle == string("e") ||
      particle == string("electron"))
    _particle = 'e';
  else if (particle == string("hl") || particle == string("h") ||
      particle == string("hole"))
    _particle = 'h';
  else
    throw ModelErrorException("Unknown particle '" + particle + "'");
}

BandProperties::~BandProperties(void)
{
  destroy(_dos_model);
}


void
BandProperties::prepare_submodels(void)
{
  //if (!get_options().has_submodel("particle_density"))
  //  get_options().add_submodel("particle_density", ModelOptions());

  //ModelOptions::submodel_iterator
  //      it(get_options().submodels_begin("particle_density"));
  //ModelOptions& opts = it->second;
  //opts["particle"] = get_option("particle","-");

  //create_submodel(_density, "particle_density", opts);


  if (!get_options().has_submodel("density_of_states"))
  {
    ModelOptions opts;
    opts.set_name("bulk");
    get_options().add_submodel("density_of_states", opts);
  }

  ModelOptions::submodel_iterator
           it(get_options().submodels_begin("density_of_states"));
  ModelOptions& dosopts = it->second;
  dosopts["particle"] = _particle;
  //if (!dosopts.find_option("level"))
  //  dosopts["level"] = get_option("band_edge", "0");

  create_submodel(_dos_model, "density_of_states", dosopts);

}



void
BandProperties::do_print_info(void)
{
  ostringstream os;
  os << "DOS: " << _dos_model->get_name();
  Messages m;
  m.info(os.str());
  m.indent();
  _dos_model->print_info();
}



double
BandProperties::get_lattice_temperature(void) const
{
  return get_driftdiffusionproperties().get_lattice_temperature();
}



bool
BandProperties::has_quantum(void) const
{
  return _dos_model->has_quantum_density();
}


void
BandProperties::use_quantum(bool use_quantum)
{
  _dos_model->use_quantum_density(use_quantum);
}


void
BandProperties::calculate(double temperature)
{
  set_temperature(temperature);
  //do_calculate();


}

double
BandProperties::get_band_edge(void) const
{
  if (_particle == 'e')
    return(*min_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
  else
    return(*max_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
}

double
BandProperties::get_effective_mass(void) const
{
  size_t i = 0;
  if (_particle == 'e')
  {
    i = distance(_dos_model->get_reference_energy().begin(),
        min_element(_dos_model->get_reference_energy().begin(),
            _dos_model->get_reference_energy().end()));
  }
  else
  {
    i = distance(_dos_model->get_reference_energy().begin(),
        max_element(_dos_model->get_reference_energy().begin(),
            _dos_model->get_reference_energy().end()));
  }
  return(_dos_model->get_effective_mass()[i]);
}

std::pair<double, double>
BandProperties::get_density_and_derivative(void) const
{
  const DriftDiffusionProperties& ddp = get_driftdiffusionproperties();

  double kT = _temperature;
  double sign = (_particle == 'h') ? 1 : -1;
  double Ef = (_particle == 'h') ?
      ddp.get_hole_electro_chemical_potential() :
      ddp.get_electron_electro_chemical_potential();
  Ef *= sign;
  double Epot = sign * ddp.get_electric_potential();

  pair<double, double> dens_der;

  if (_dos_model->has_quantum_density() && ddp.has_solution()) // && use_predictor())
  {
    double Ef_old = (_particle == 'h') ? ddp.get_old_fermih() : -ddp.get_old_fermie();
    double Epot_old = sign * ddp.get_old_phi();

    // get the old quantum density
    pair<double, double> dens_der_old(_dos_model->get_occupied_density_and_derivative(
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature()));

    // now get the old classical density
    _dos_model->use_quantum_density(false);

    pair<double, double> dens_der_old_cl(_dos_model->get_occupied_density_and_derivative(
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature()));

    // now get the new classical density and derivative
    dens_der = _dos_model->get_occupied_density_and_derivative(
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());

    double fac = dens_der_old.first / dens_der_old_cl.first;
    dens_der.first *= fac;
    dens_der.second *= -sign * fac;

    _dos_model->use_quantum_density(true);
  }
  else
  {
    dens_der = _dos_model->get_occupied_density_and_derivative(
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());
    dens_der.second *= -sign;
  }

  return dens_der;
}



std::pair<double, double>
BandProperties::get_density_and_derivative(double Ef, double Epot) const
{
    const DriftDiffusionProperties& ddp = get_driftdiffusionproperties();

  double kT = _temperature;
  double sign = (_particle == 'h') ? 1 : -1;
  //double Ef = (_particle == 'h') ?
  //    ddp.get_hole_electro_chemical_potential() :
  //    ddp.get_electron_electro_chemical_potential();
  //Ef *= sign;
  //double Epot = sign * ddp.get_electric_potential();
  Epot *= sign;
  Ef *= sign;

  pair<double, double> dens_der;

  if (_dos_model->has_quantum_density() && ddp.has_solution()) // && use_predictor())
  {
    double Ef_old = (_particle == 'h') ? ddp.get_old_fermih() : -ddp.get_old_fermie();
    double Epot_old = sign * ddp.get_old_phi();

    // get the old quantum density
    pair<double, double> dens_der_old(_dos_model->get_occupied_density_and_derivative(
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature()));

    // now get the old classical density
    _dos_model->use_quantum_density(false);

    pair<double, double> dens_der_old_cl(_dos_model->get_occupied_density_and_derivative(
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature()));

    // now get the new classical density and derivative
    dens_der = _dos_model->get_occupied_density_and_derivative(
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());

    double fac = dens_der_old.first / dens_der_old_cl.first;
    dens_der.first *= fac;
    dens_der.second *= -sign * fac;

    _dos_model->use_quantum_density(true);
  }
  else
  {
    dens_der = _dos_model->get_occupied_density_and_derivative(
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());
    dens_der.second *= -sign;
  }

  return dens_der;



  /*

  const DriftDiffusionProperties& ddp = get_driftdiffusionproperties();
  _density->set_element_and_point(ddp.get_element(), ddp.get_coordinates());

  double edge = get_band_edge();
  double kT = _temperature;
  double sign = (_particle == 'h') ? -1 : 1;

  double dens = 0;
  double dens_der = 0;

  if (_density->is_quantum_density() && ddp.has_solution()) // && use_predictor())
  {
    double Ef_old = (_particle == 'h') ? ddp.get_old_fermih() : -ddp.get_old_fermie();

    // set the OLD potentials
    _density->set_classical_parameters(get_effective_mass(),
        sign * (edge - ddp.get_old_phi()), Ef_old, kT);

    double qdens = _density->get_particle_density();

    // now get the old classical density
    _density->use_quantum_density(false);
    double old_dens = _density->get_particle_density();

    // now get the new classical density and derivative
    _density->set_classical_parameters(get_effective_mass(),
        sign * edge + Epot, Ef, kT);
    dens = _density->get_particle_density();
    dens_der = sign * _density->get_particle_density_derivative();

    double fac = qdens / old_dens;
    dens *= fac;
    dens_der *= fac;

    _density->use_quantum_density(true);
  }
  else
  {
    / *
    _density->set_classical_parameters(get_effective_mass(),
    sign * edge + Epot, Ef, kT);
    dens = _density->get_particle_density();
    dens_der = sign * _density->get_particle_density_derivative();
    * /


    pair<double, double> dens_der(_dos_model->get_occupied_density_and_derivative(
    Ef, Epot, kT, ddp.get_element(), ddp.get_coordinates(), ddp.get_lattice_temperature() ) );
    dens_der.second *= sign;
    return dens_der;

  }

  return std::make_pair(dens, dens_der);
  */
}

double
BandProperties::get_gamma(void) const
{
  return 1; //_density->get_gamma();
};
