// $Id$

#include "MyPoissonModel.h"
#include "MyChargeDensityModel.h"


using namespace std;


TiberModelObject*
MyPoissonModel::_create(const ModelOptions& options)
{
  return new MyPoissonModel(options);
}


void
MyPoissonModel::_destroy(TiberModelObject* p)
{
  delete p;
}


MyPoissonModel*
MyPoissonModel::create(const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  PhysicalModelInterface* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModelInterface::create(_create, _destroy, options);
  else
  {
    type = "bulk_" + type;
    pm = PhysicalModelInterface::create(type, options);
  }

  return dynamic_cast<MyPoissonModel*>(pm);
}



void
MyPoissonModel::do_init(void)
{
  SubmodelIterator it = submodels_begin("charge_density");
  if (it != submodels_end("charge_density"))
  {
    _charge_density = dynamic_cast<ChargeDensityModel*>(it->second);
  }
}



void
MyPoissonModel::do_calculate(void)
{
  _permittivity(0,0) = _permittivity(1,1) = _permittivity(2,2) = 1.0;

  _polarization.zero();

  _charge = 0.0;
  if (_charge_density != NULL)
    _charge = _charge_density->get_charge_density(get_element(), get_point());

}
