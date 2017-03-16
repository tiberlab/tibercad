// $Id$

#include "DensityOfStates.h"
#include "ModelErrorException.h"

#include <string>

using namespace std;

DensityOfStates::DensityOfStates(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _fixed_DOS(false),
  _reference_energy({0.0}),
  _effective_mass({1.0}),
  _particle(' '),
  _spin(0.5),
  _use_quantum(false),
  _is_quantum(false)
{
  string particle = get_option("particle", "-");
  if (particle == string("el") || particle == string("e") ||
      particle == string("electron"))
    _particle = 'e';
  else if (particle == string("hl") || particle == string("h") ||
      particle == string("hole"))
    _particle = 'h';

  _spin = get_option("spin", _spin);
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
      PhysicalModelInterface::create("density_of_states_" + name, NULL, options));
  if (dos_ptr == NULL)
    throw ModelErrorException("Unknown density of states type: " + name);

  return dos_ptr;
}


