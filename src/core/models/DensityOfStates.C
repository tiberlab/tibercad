// $Id$

#include "DensityOfStates.h"
#include "ExponentialDOS.h"
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

  string name(options.get_key());
  name = options.get_option("type", name);
  if (name.empty())
    name = "delta";

  if (name == "delta")
    dos_ptr = NULL;
  else if (name == "exponential")
    dos_ptr = ExponentialDOS::create(options);
  else
    throw ModelErrorException("Unknown density of states type: " + name);

  return dos_ptr;
}


