// $Id$

#include "CarrierProperties.h"
#include "ParticleDensity.h"
#include "DensityOfStates.h"
#include "MobilityModelInterface.h"
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
    _background_conductivity(1e-12),
    _conductivity_model(STD),
    _conductivity(0.0),
    _conductivity_derivative_qFermi(0.0),
    _conductivity_derivative_gradqFermi(0.0),
    _conductivity_derivative_gradPotential(0.0),
    _dos_factor(pow(2.0 * M_PI *
        Constants::me / (Constants::h * Constants::h) *
        Constants::e, 1.5) / 1e6),
    _dos_model(nullptr),
    _mobility_model(nullptr)
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
  destroy(_mobility_model);
}


void
CarrierProperties::do_init(void)
{
  _carrier_id = this->get_driftdiffusionproperties().get_carrier_id(_particle_name);

  _background_conductivity = get_option("background_conductivity", _background_conductivity);

  string conductivity_model("standard");
  conductivity_model = get_option("conductivity_model", conductivity_model);
  if (conductivity_model == "standard")
    _conductivity_model = STD;
  else if (conductivity_model == "generalized")
    _conductivity_model = GEN;
  else if (conductivity_model == "constant")
  {
    _conductivity_model = CONST;
    // TODO read from database, or implement conductivity models
    _background_conductivity = 1.0;
    get_parameter("conductivity", _background_conductivity);
  }
  else
    throw InitFailedException("No such conductivity model: " + conductivity_model);
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

  if (!get_options().has_submodel("mobility"))
  {
    ModelOptions opts;
    opts.set_name("constant");
    get_options().add_submodel("mobility", opts);
  }
  it = get_options().submodels_begin("mobility");

  create_submodel(_mobility_model, "mobility", it->second);
  _mobility_model->set_carrier_type(_particle);
  _mobility_model->set_carrier(_carrier_id);

  if (++it != get_options().submodels_end("mobility"))
    throw InitFailedException("Carrier " + this->get_particle_name() +
        " : cannot have more than one mobility model per region");

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
  _mobility_model->print_info();
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

double
CarrierProperties::get_mobility(void) const
{
  double mobility = 1.0;
  if (_mobility_model != nullptr)
  {
    mobility = _mobility_model->get_mobility();
  }

  return(mobility);
}


std::pair<double, double>
CarrierProperties::get_density_and_derivative(double Ef, double Epot)
{
  const DriftDiffusionProperties& ddp = get_driftdiffusionproperties();

  double kT = _temperature;
  double sign = (_charge > 0.0) ? 1 : -1;
  Epot *= _charge;
  Ef *= sign;


  vector<double> dens_der(2);
  if (_conductivity_model == GEN)
    dens_der.resize(3);

  if (_dos_model->has_quantum_density() && ddp.has_solution()) // && use_predictor())
  {
    double Ef_old = ddp.get_old_q_fermi_potential(_carrier_id);
    Ef_old *= sign;
    double Epot_old = _charge * ddp.get_old_phi();

    // get the old quantum density
    auto dens_der_old(dens_der);
    _dos_model->get_occupied_density_and_derivative(dens_der_old,
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());

    // now get the old classical density
    _dos_model->use_quantum_density(false);

    auto dens_der_old_cl(dens_der);
    _dos_model->get_occupied_density_and_derivative(dens_der_old_cl,
          Ef_old, Epot_old, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());

    // now get the new classical density and derivative
    _dos_model->get_occupied_density_and_derivative(dens_der,
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());

    double fac = dens_der_old[0] / dens_der_old_cl[0];
    dens_der[0] *= fac;
    dens_der[1] *= -sign * fac;
    if (_conductivity_model == GEN)
      dens_der[2] *= sign * sign * fac;


    _dos_model->use_quantum_density(true);
  }
  else
  {
    _dos_model->get_occupied_density_and_derivative(dens_der,
          Ef, Epot, kT, ddp.get_element(),
          ddp.get_coordinates(), ddp.get_lattice_temperature());
    dens_der[1] *= -sign;
    if (_conductivity_model == GEN)
      dens_der[2] *= sign * sign;
  }


  double mobility = _mobility_model->get_mobility();
  switch (_conductivity_model)
  {
    case STD:
      _conductivity = mobility * dens_der[0];
      _conductivity_derivative_qFermi = mobility * dens_der[1];

      _mobility_model->get_derivative_grad_fermi(_conductivity_derivative_gradqFermi);
      _conductivity_derivative_gradqFermi *= dens_der[0];

      _mobility_model->get_derivative_grad_potential(_conductivity_derivative_gradPotential);
      _conductivity_derivative_gradPotential *= dens_der[0];
      break;

    case GEN:

      mobility *= kT;
      _conductivity = -sign * mobility * dens_der[1];
      _conductivity_derivative_qFermi = -sign * mobility * dens_der[2];

      _mobility_model->get_derivative_grad_fermi(_conductivity_derivative_gradqFermi);
      _conductivity_derivative_gradqFermi *= -sign * dens_der[1];

      _mobility_model->get_derivative_grad_potential(_conductivity_derivative_gradPotential);
      _conductivity_derivative_gradPotential *= -sign * dens_der[1];
      break;

    default:
      break;
  }

  _conductivity += _background_conductivity / Constants::e;

  return(make_pair(dens_der[0], dens_der[1]));
}
