// $Id$

#include "Multiscale.h"
#include "ModelOptions.h"
#include "Messages.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "InitFailedException.h"

#include <algorithm>

using namespace std;

Multiscale::Multiscale(const ModelOptions& options)
: _options(options)
{
  init();
}



void
Multiscale::init(void)
{
  string method(_options.get_name());
  method = _options.get_option("type", method);

  if (method == "bridge")
    _method = BRIDGE;
  else if (method == "overlap")
    _method = OVERLAP;

  _macromodel = SimulationInterface::find_simulation(_options.get_option("macromodel", ""));
  _micromodel = SimulationInterface::find_simulation(_options.get_option("micromodel", ""));
  if ((_macromodel == NULL) || (_micromodel == NULL))
    throw InitFailedException("Multiscale definition needs both macro- and microscale models");

  if (!(_macromodel->has_environment() && _micromodel->has_environment()))
    throw InitFailedException("Multiscale definition needs both macro- and microscale models "
        "with associated mesh regions");

  ostringstream os;
  os << "Using multiscale method: ";
  if (_method == BRIDGE)
  {
    os << "bridge\n";
    // find domain restriction for macroscale model
    _options.get_option("restrict_variables", _restricted_variables);

    set<ID> macro_ids, micro_ids;
    _macromodel->get_region_ids(macro_ids);
    _micromodel->get_region_ids(micro_ids);

    vector<ID> diff(macro_ids.size());
    vector<ID>::iterator last = set_difference(macro_ids.begin(), macro_ids.end(),
        micro_ids.begin(), micro_ids.end(),
        diff.begin());

    set<string> regnames;
    for (int i = 0; i < last - diff.begin(); ++i)
    {
      _restricted_macro_domain.insert(diff[i]);
      regnames.insert(_macromodel->get_environment().get_device().get_region_name(diff[i]));
    }

    os << "  macromodel restricted to regions:";
    for (set<string>::iterator it(regnames.begin()); it != regnames.end(); ++it)
      os << " " << *it;
  }
  else
    os << "overlap";

  Messages::newline();
  Messages::info(os.str());

}
