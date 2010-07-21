// $Id$

#include "MyPoissonModel.h"


using namespace std;



MyPoissonModel*
MyPoissonModel::create(const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);
  type = "bulk_" + type;

  return dynamic_cast<MyPoissonModel*>(
      PhysicalModelInterface::create(type, options));
}



void
MyPoissonModel::do_init(void)
{
  SubmodelIterator it = submodels_begin("charge_density");
  if (it != submodels_end("charge_density"))
    _charge_density = dynamic_cast<ChargeDensityModel*>(it->second);
}



void
MyPoissonModel::do_calculate(void)
{
  _permittivity(0,0) = _permittivity(1,1) = _permittivity(2,2) = 1.0;

  _polarization.zero();

  _charge = 0.0;
  if (_charge_density != NULL)
    charge = _charge_density->get_charge_density(get_element(), get_point());

}
