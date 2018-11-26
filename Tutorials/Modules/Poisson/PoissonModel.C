// $Id$

#include "PoissonModel.h"
#include "ChargeDensityModel.h"
#include "PolarizationModel.h"
#include "PermittivityModel.h"

using namespace std;


TiberModelObject*
PoissonModel::_create(const ModelOptions& options, const void*)
{
  return new PoissonModel(options);
}


void
PoissonModel::_destroy(TiberModelObject* p)
{
  delete p;
}


PoissonModel*
PoissonModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  PoissonModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModelInterface::create<PoissonModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    type = "bulk_" + type;
    pm = PhysicalModelInterface::create<PoissonModel>(type, mat, options);
  }

  return(pm);
}



void
PoissonModel::do_init(void)
{
  _permittivity = _permittivity_model->get_permittivity();
}

void
PoissonModel::do_calculate(void)
{
 

  _charge_density = 0.0;
  if (_charge_density_model != NULL)
    _charge_density = _charge_density_model->get_charge_density(get_element(), get_point());

  //Get total polarization 
  _polarization = 0;
  for (ID n = 0 ; n <_pm.size() ; n++)
  {
    _pm[n]->calculate(_elem,_point);
    _polarization += _pm[n]->get_polarization();
  }

  

}

void
PoissonModel::prepare_submodels(void)
{
  
  ModelOptions opts;
  opts.set_option("type", "constant");
  create_submodel(_permittivity_model, "permittivity", opts);
  
  // alternative way to create internal submodels:
  //
  // PermittivityModel* mod = PhysicalModelInterface::create("permittivity", opts);
  // add_submodel("permittivity", mod)

  create_submodel(_charge_density_model, "charge_density");
  create_submodels(_pm, "polarization");

  // NOTE: all submodels are initialized automatically before calling do_init()
}
