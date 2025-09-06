#include "DegradationH2O.h"
#include "SimulationInterface.h"

#include "TiberModule.h"


DegradationH2O::DegradationH2O(const ModelOptions& options)
  : DegradationModel(options)
{
}

DegradationH2O*
DegradationH2O::create(const ModelOptions& options)
{
  return new DegradationH2O(options);
}

void
DegradationH2O::do_init(void)
{

  //_RH_ref = get_option("reference_humidity", _RH_ref);

  //_exponent = get_option("exponent", _exponent);

  std::string water_ingress_sim = get_option("relative_humidity", "");
  _humidity_model = SimulationInterface::find_solution_provider(water_ingress_sim, "RelativeHumidity");
}



