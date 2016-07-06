// $Id$

#include "MEBandProperties.h"
#include "ParticleDensity.h"
#include "DensityOfStates.h"
#include "MasterEquationsProperties.h"
#include "ModelOptions.h"
#include "Messages.h"
#include "Database.h"


#include <sstream>
#include <algorithm>

using namespace std;



MEBandProperties::MEBandProperties(const ModelOptions& options) :
    MasterEquationsModelInterface(options),
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

MEBandProperties::~MEBandProperties(void)
{
  destroy(_dos_model);
}


//void
//MEBandProperties::prepare_submodels(void)
//{

 // if (!get_options().has_submodel("density_of_states"))
  //{
    //ModelOptions opts;
    //opts.set_name("me_bulk");
    //get_options().add_submodel("density_of_states", opts);
  //}
//
  //ModelOptions::submodel_iterator
           //it(get_options().submodels_begin("density_of_states"));
  //ModelOptions& dosopts = it->second;
  //dosopts["particle"] = _particle;
  //if (!dosopts.find_option("level"))
  //  dosopts["level"] = get_option("band_edge", "0");

  //create_submodel(_dos_model, "density_of_states", dosopts);

//}



void
MEBandProperties::do_print_info(void)
{
  ostringstream os;
  os << "DOS: " << _dos_model->get_name();
  Messages m;
  m.info(os.str());
  m.indent();
  _dos_model->print_info();
}



double
MEBandProperties::get_lattice_temperature(void) const
{
  return get_masterequationsproperties().get_lattice_temperature(); // ??
}



void
MEBandProperties::calculate(double temperature)
{
  set_temperature(temperature);
  //do_calculate();


}

double
MEBandProperties::get_band_edge(void) const
{
  if (_particle == 'e')
    return(*min_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
  else
    return(*max_element(_dos_model->get_reference_energy().begin(),
        _dos_model->get_reference_energy().end()));
}

double
MEBandProperties::get_effective_mass(void) const
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

