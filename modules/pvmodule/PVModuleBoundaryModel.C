// $Id$

#include "PVModuleBoundaryModel.h"
#include "Messages.h"

using namespace std;


PVModuleBoundaryModel::PVModuleBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


TiberModelObject*
PVModuleBoundaryModel::_create(const ModelOptions& options, const void*)
{
  return new PVModuleBoundaryModel(options);
}


void
PVModuleBoundaryModel::_destroy(TiberModelObject* p)
{
  delete p;
}


PVModuleBoundaryModel*
PVModuleBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  PVModuleBoundaryModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<PVModuleBoundaryModel>(_create, _destroy, boundary, options);
  else
  {
    // there is no such model, at the moment
    type = "contact_" + type;
    pm = PhysicalModel::create<PVModuleBoundaryModel>(type, boundary, options);
  }

  return(pm);
}



void
PVModuleBoundaryModel::do_init(void)
{
  string type = get_option("type", "ground");
  if (type == "ground") _contact_type = GND;
  if (type == "source") _contact_type = SRC;

  string layer = get_option("layer", "bottom");
  if (layer == "bottom") _contact_layer = BOTTOM;
  if (layer == "top") _contact_layer = TOP;
  if (layer == "both") _contact_layer = BOTH;
}

