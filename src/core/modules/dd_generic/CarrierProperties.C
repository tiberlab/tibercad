// $Id$

#include "CarrierProperties.h"
#include "ParticleDensity.h"
#include "DensityOfStates.h"
#include "DriftDiffusionProperties.h"
#include "ModelOptions.h"
#include "Messages.h"
#include "Database.h"


#include <sstream>
#include <algorithm>

using namespace std;



CarrierProperties::CarrierProperties(const ModelOptions& options) :
    DriftDiffusionModelInterface(options),
    _particle('e'),
    _particle_name(""),
    _carrier_id(DriftDiffusionProperties::unknown_carrier_id),
    _charge(-1),
    _spin(0.5),
    _is_dopant(false),
    _dos_factor(pow(2.0 * M_PI *
        Constants::me / (Constants::h * Constants::h) *
        Constants::e, 1.5) / 1e6)
{
  _particle_name = get_options().get_name();
  _particle_name = get_option("name", _particle_name);

  if (_particle_name.empty())
    throw ModelErrorException("Carrier name MUST be provided");

  if (_particle_name == string("el") ||
      _particle_name == string("e") ||
      _particle_name == string("electron") ||
      _particle_name == string("electrons") ||
      _particle_name == string("n") ||
      _particle_name == string("negative"))
  {
    _charge = -1;
  }
  else if (_particle_name == string("hl") ||
      _particle_name == string("h") ||
      _particle_name == string("hole") ||
      _particle_name == string("holes") ||
      _particle_name == string("p") ||
      _particle_name == string("positive"))
  {
    _charge = 1;
  }
  else if (_particle_name == string("ex") ||
      _particle_name == string("x") ||
      _particle_name == string("exciton"))
  {
    _charge = 0;
    _spin = 0;
  }

  _charge = get_option("charge", _charge);
  _spin = get_option("spin", _spin);

  _particle = (_charge <= 0.0) ? 'e' : 'h';

  _is_dopant = get_option("dopant", _is_dopant);

}

CarrierProperties::~CarrierProperties(void)
{
  destroy(_dos_model);
}


void
CarrierProperties::do_init(void)
{
  _carrier_id = this->get_driftdiffusionproperties().get_carrier_id(_particle_name);
}

void
CarrierProperties::do_reinit(void)
{
}

void
CarrierProperties::prepare_submodels(void)
{
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
  dosopts["charge"] = _charge;
  dosopts["spin"] = _spin;

  //if (!dosopts.find_option("level"))
  //  dosopts["level"] = get_option("band_edge", "0");


  create_submodel(_dos_model, "density_of_states", dosopts);

  //for an exciton reference_energy() will store the band edges of the corresponding e and h
  /*
  if (_is_exciton)
  {
    _dos_model->reference_energy().resize(0);
    _dos_model->reference_energy().push_back(0.0);  // reserve an empty place for the dos model to set the ref energy as Eg - R (R being the binding energy of the exciton)

    for (size_t i = 0; i<_exciton_carriers.size(); i++)
      _dos_model->reference_energy().push_back(  ( get_driftdiffusionproperties().get_carrier_properties(
                                                     _exciton_carriers[i]) )->get_band_edge()  );
  }
  */
}



void
CarrierProperties::do_print_info(void)
{
  ostringstream os;
  os << "name = " << _particle_name
     << ", charge = " << _charge
     << ", spin = " << _spin << "\n";
  os << "DOS: " << _dos_model->get_name();
  Messages m;
  m.info(os.str());
  m.indent();
  _dos_model->print_info();
}


double
CarrierProperties::get_maximum_density(void) const
{
  return(_dos_model->get_total_state_density());
}


double
CarrierProperties::get_lattice_temperature(void) const
{
  return get_driftdiffusionproperties().get_lattice_temperature();
}



bool
CarrierProperties::has_quantum(void) const
{
  return(const_cast<const DensityOfStates*>(_dos_model)->is_quantum_density());
}


void
CarrierProperties::use_quantum(bool use_quantum)
{
  //_dos_model->use_quantum_density(use_quantum);
}


void
CarrierProperties::calculate(double temperature)
{
  set_temperature(temperature);
  //do_calculate();
}

double
CarrierProperties::get_band_edge(void) const
{
  if (_charge <= 0)
    return(*min_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
  else
    return(*max_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
}

double
CarrierProperties::get_effective_mass(void) const
{
  size_t i = 0;
  if (_charge <= 0)
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


double
CarrierProperties::get_thermoelectric_power(void) const
{
  return _dos_model->get_thermoelectric_power();
}


std::pair<double, double>
CarrierProperties::get_density_and_derivative(double Ef, double Epot) const
{
  const DriftDiffusionProperties& ddp = get_driftdiffusionproperties();

  double kT = _temperature;
  double sign = (_charge > 0.0) ? 1 : -1;
  Epot *= _charge;
  Ef *= sign;


  pair<double, double> dens_der;

  if (_dos_model->has_quantum_density() && ddp.has_solution()) // && use_predictor())
  {
    double Ef_old = ddp.get_old_q_fermi_potential(_carrier_id);
    Ef_old *= sign;
    double Epot_old = _charge * ddp.get_old_phi();

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
