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
  SubmodelIterator it = submodels_begin("charge_density");
  if (it != submodels_end("charge_density"))
    _charge_density = dynamic_cast<MyChargeDensityModel*>(it->second);
  
  it = submodels_begin("polarization");
  const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("polarization"));
  for ( ; it != it_end ; ++it)
    _pm.push_back(dynamic_cast<PolarizationModel*> ((*it).second));

  it = submodels_begin("permittivity");
  if (it != submodels_end("permittivity"))
    _permittivity_model = dynamic_cast<PermittivityModel*>(it->second);

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
  
  //Thermal Conductivity Default
  if (!get_options().has_submodel("permittivity"))
  {
    ModelOptions opts;
    opts.set_option("type","constant");
    get_options().add_submodel("permittivity",opts);
  }
  
}
