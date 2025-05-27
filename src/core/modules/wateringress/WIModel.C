// $Id$

#include "WIModel.h"

using namespace std;


TiberModelObject*
WIModel::_create(const ModelOptions& options, const void*)
{
  return new WIModel(options);
}


void
WIModel::_destroy(TiberModelObject* p)
{
  delete p;
}


WIModel*
WIModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  WIModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<WIModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    type = "bulk_" + type;
    pm = PhysicalModel::create<WIModel>(type, mat, options);
  }

  return(pm);
}



void
WIModel::do_init(void)
{
  // we read it in g/m^3/Pa
  _solubility = get_option("solubility", _solubility);
  _solubility / 1e6;

  // we read it in M^2/S
  _diffusivity = get_option("diffusivity", _diffusivity);
  _diffusivity * 1e4;
}


void
WIModel::calculate(const Elem* elem, const Point& point)
{
 
}

void
WIModel::prepare_submodels(void)
{
}
