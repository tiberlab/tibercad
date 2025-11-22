// $Id$

#include "tibercad/physics/misc/DensityOfStates.h"
#include "tibercad/base/ModelErrorException.h"
#include "tibercad/profiles/ExternalProfile.h"

#include <string>

using namespace std;

DensityOfStates::DensityOfStates(const ModelOptions& options) :
  PhysicalModel(options),
  _fixed_DOS(false),
  _reference_energy({0.0}),
  _effective_mass({1.0}),
  _effective_dos(1e19),
  _particle(' '),
  _spin(0.5),
  _use_quantum(false),
  _is_quantum(false),
  _th_el_power(0.0),
  _total_density(1e23),
  _profile(nullptr)
{
  string particle = get_option("particle", "-");
  if (particle == string("el") || particle == string("e") ||
      particle == string("electron"))
    _particle = 'e';
  else if (particle == string("hl") || particle == string("h") ||
      particle == string("hole"))
    _particle = 'h';

  _spin = get_option("spin", _spin);

  _effective_dos = get_option("N0", _effective_dos);

  if (get_options().has_submodel("profile"))
  {
    _profile = ExternalProfile::create(
        get_options().submodels_begin("profile")->second);
  }
}


DensityOfStates*
DensityOfStates::create(const ModelOptions& options)
{
  DensityOfStates* dos_ptr;

  string name(options.get_name());
  name = options.get_option("type", name);
  if (name.empty())
    name = "delta";

  dos_ptr = dynamic_cast<DensityOfStates*>(
      PhysicalModel::create("density_of_states_" + name, NULL, options));
  if (dos_ptr == NULL)
    throw ModelErrorException("Unknown density of states type: " + name);

  return dos_ptr;
}

double
DensityOfStates::get_effective_dos(const Elem* elem, const Point& p) const
{
  double dos = _effective_dos;
  if (_profile != nullptr)
    dos *= _profile->get_data(elem, p);

  return(dos);
}

