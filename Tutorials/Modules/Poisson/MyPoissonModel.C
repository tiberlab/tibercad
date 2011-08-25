// $Id$

#include "MyPoissonModel.h"
#include "MyChargeDensityModel.h"
#include "PolarizationModel.h"
#include "PermittivityModel.h"

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
  _permittivity = _permittivity_model->get_permittivity();
}

void
MyPoissonModel::do_calculate(void)
{
 

  _charge = 0.0;
  if (_charge_density != NULL)
    _charge = _charge_density->get_charge_density(get_element(), get_point());

  //Get total polarization 
  _polarization = 0;
  for (ID n = 0 ; n <_pm.size() ; n++)
  {
    _pm[n]->calculate(_elem,_point);
    _polarization += _pm[n]->get_polarization();
  }

  

}

void
MyPoissonModel::prepare_submodels(void)
{
  
  ModelOptions opts;
  opts.set_option("type", "constant");
  create_submodel(_permittivity_model, "permittivity", opts);
  
  // alternative way to create internal submodels:
  //
  // PermittivityModel* mod = PhysicalModelInterface::create("permittivity", opts);
  // add_submodel("permittivity", mod)

  create_submodel(_charge_density, "charge_density");
  create_submodels(_pm, "polarization");

  // NOTE: all submodels are initialized automatically before calling do_init()
}
