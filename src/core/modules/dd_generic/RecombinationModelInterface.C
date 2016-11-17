// $Id$
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "Material.h"

void
RecombinationModelInterface::do_init(void)
{

  get_option("carriers", _carriers);

  if (_carriers.size() < 2)
    throw InitFailedException("Recombination model '" + get_default_name() + "' requires to provide two carriers");

  _carriers.resize(2);

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  for ( auto it : _carriers)
  {
    if (dd.get_carrier_properties(it) == nullptr)
      throw InitFailedException("Recombination '" + get_default_name() + "': carrier '" + (it) + "' not found'");
  }

  if (_carriers[0] == _carriers[1])
    throw InitFailedException("Recombination '" + get_default_name() + "': carrier names must be different");
}