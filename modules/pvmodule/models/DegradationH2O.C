#include "DegradationH2O.h"
#include "tibercad/module/SimulationInterface.h"

#include "tibercad/module/TiberModule.h"


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

  if (get_options().has_submodel("photocurrent"))
  {
    ModelOptions& opts = get_options().submodels_begin("photocurrent")->second;
    _RH_ref_ph = opts.get_option("reference_humidity", _RH_ref_ph);
    _exponent_ph = opts.get_option("exponent", _exponent_ph);
  }

  if (get_options().has_submodel("series_resistance"))
  {
    ModelOptions& opts = get_options().submodels_begin("series_resistance")->second;
    _RH_ref_rs = opts.get_option("reference_humidity", _RH_ref_rs);
    _exponent_rs = opts.get_option("exponent", _exponent_rs);
  }
 
  std::string water_ingress_sim = get_option("relative_humidity", "");
  _humidity_model = SimulationInterface::find_solution_provider(water_ingress_sim, "RelativeHumidity");
}

void
DegradationH2O::do_degrade_params(const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  DegradationModel::Parameters& params) const
{
  if (_humidity_model.first != nullptr)
  {
    double humidity = 0;

    _humidity_model.first->get_solution(elem, _humidity_model.second, humidity, p, false);

    // factor for photocurrent
    double ph_fac = 1.0 + std::pow(humidity / _RH_ref_ph, _exponent_ph);

    params.double_params[2] /= ph_fac;

    // factor for series resistance
    double rs_fac = (1.0 + std::pow(humidity / _RH_ref_rs, _exponent_rs));

    params.double_params[0] *= rs_fac;
  }
}
