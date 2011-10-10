// $Id$

#include "DensityOfStates.h"
#include "ModelErrorException.h"

#include <string>

using namespace std;

DensityOfStates::DensityOfStates(const ModelOptions& options) :
    PhysicalModelInterface(options)
{
}


DensityOfStates*
DensityOfStates::create(const ModelOptions& options)
{
  DensityOfStates* dos_ptr;

  string name(options.get_name());
  name = options.get_option("type", name);
  if (name.empty())
    name = "delta";

  if (name == "delta")
    dos_ptr = NULL;
  else
  {
    dos_ptr = dynamic_cast<DensityOfStates*>(
        PhysicalModelInterface::create("density_of_states_" + name, NULL, options));
    if (dos_ptr == NULL)
      throw ModelErrorException("Unknown density of states type: " + name);
  }

  return dos_ptr;
}


