// $Id$
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "Material.h"

using namespace std;


void
RecombinationModelInterface::do_init(void)
{

  _plot_name = get_option("plot_name", "");
  get_option("carriers", _carriers);
  get_option("exponents", _exponents);

  _is_radiative = get_option("radiative", _is_radiative);

  // to be sure
  _carrier_ids.resize(0);

  // we are paranoid and check for multiple names
  set<ID> used_ids;

  for ( auto it : _carriers)
  {
    ID id = get_bulk_driftdiffusionproperties().get_carrier_id(it);
    if (id == DriftDiffusionProperties::unknown_carrier_id)
      throw InitFailedException("Recombination '" + get_default_name() +
          "': carrier '" + (it) + "' not found in module");

    //const CarrierProperties* cp = get_bulk_driftdiffusionproperties().get_carrier_properties(id);
    const CarrierProperties* cp = get_driftdiffusionproperties().get_carrier_properties(id);
    //if ( cp == nullptr)
    //  throw InitFailedException("Recombination '" + get_default_name() +
    //      "': carrier '" + (it) + "' not defined in the same regions as the recombination model");

    if (used_ids.count(id))
      throw InitFailedException("Recombination '" + get_default_name() +
          "': carrier names must be different");

    _carrier_ids.push_back(id);

    used_ids.insert(id);
  }

  if (_exponents.size() == 0)
  {
    _exponents.resize(_carriers.size(), 1);
  }

  if (_exponents.size() != _carriers.size())
  {
    throw InitFailedException("Recombination '" + get_default_name() +
        "': number of exponents inconsistent with number of carriers");
  }

  // for the case of electrons and holes we put them in a standard order
  // electron, hole
  if (_carriers.size() == 2)
  {
    if ((_carriers[1] == "electron") && (_carriers[0] == "hole"))
    {
      vector<ID> new_order(2);
      new_order[0] = _carrier_ids[1];
      new_order[1] = _carrier_ids[0];
      this->reorder_carriers(new_order);
    }
  }
}



void
RecombinationModelInterface::reorder_carriers(const std::vector<ID>& new_order)
{
  vector<string> old_order(_carriers);

  assert(new_order.size() == _carriers.size());

  for (unsigned int i = 0; i < new_order.size(); ++i)
  {
    for (unsigned int j = 0; j < _carriers.size(); ++j)
    {
      if (_carrier_ids[j] == new_order[i])
      {
        _carriers[i] = old_order[j];
      }
    }
  }

  _carrier_ids = new_order;
}

void
RecombinationModelInterface::reorder_ids(const std::vector<std::string>& new_order)
{
  vector<ID> old_order(_carrier_ids);

  assert(new_order.size() == _carrier_ids.size());

  for (unsigned int i = 0; i < new_order.size(); ++i)
  {
    for (unsigned int j = 0; j < _carrier_ids.size(); ++j)
    {
      if (_carriers[j] == new_order[i])
      {
        _carrier_ids[i] = old_order[j];
      }
    }
  }

  _carriers = new_order;
}

void
RecombinationModelInterface::get_net_rate_and_derivatives(
    std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  R.resize(this->get_driftdiffusionproperties().n_known_carriers(), 0.0);
  dPotentials.resize(this->get_driftdiffusionproperties().n_known_carriers(),
      vector<double>(this->get_driftdiffusionproperties().n_known_carriers() + 1, 0.0) );
  calculate_rate_and_derivatives(R, dPotentials);
}

void
RecombinationModelInterface::calculate_rate_and_derivatives(
    std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

}
