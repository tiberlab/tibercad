// $Id$
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "Material.h"

void
RecombinationModelInterface::do_init(void)
{

  get_option("carriers", _carriers);

  // to be sure
  _carrier_ids.resize(0);

  // we are paranoid and check for multiple names
  std::set<ID> used_ids;

  for ( auto it : _carriers)
  {
    ID id = get_driftdiffusionproperties().get_carrier_id(it);
    if (id == DriftDiffusionProperties::unknown_carrier_id)
      throw InitFailedException("Recombination '" + get_default_name() +
          "': carrier '" + (it) + "' not found in module");

    if (used_ids.count(id))
      throw InitFailedException("Recombination '" + get_default_name() +
          "': carrier names must be different");

    _carrier_ids.push_back(id);
    used_ids.insert(id);

  }
}


double
RecombinationModelInterface::calculate_rate_and_derivatives(
    std::vector<double>& dPotentials)
{
  dPotentials.resize(0);
  dPotentials.resize(this->get_driftdiffusionproperties().get_known_carriers().size(), 0.0);
  return(0.0);
}
