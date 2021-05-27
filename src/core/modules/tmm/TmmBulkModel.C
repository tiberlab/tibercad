// $Id: TmmBulkModel.C 4391 2017-04-07 11:16:58Z pecchia $

#include "TmmBulkModel.h"
//#include "PermittivityModel.h"

using std::string;
using namespace libMesh;


TiberModelObject*
TmmBulkModel::_create(const ModelOptions& options, const void*)
{
  return new TmmBulkModel(options);
}


inline
TmmBulkModel::TmmBulkModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}



void
TmmBulkModel::_destroy(TiberModelObject* p)
{
  delete p;
}


TmmBulkModel*
TmmBulkModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  TmmBulkModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<TmmBulkModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    //type = "bulk_" + type;
    //pm = PhysicalModel::create<TmmBulkModel>(type, mat, options);
  }

  return(pm);
}



libMesh::Complex
TmmBulkModel::get_permittivity(double lambda) const
{
  Complex eps = 1;

  return(eps);
}


void
TmmBulkModel::do_init(void)
{
}

void
TmmBulkModel::prepare_submodels(void)
{
  // Maybe it would be more elegant to extend the existing
  // permittivity model implementation
  //ModelOptions opts;
  //opts.set_option("type", "constant");
  //create_submodel(_permittivity_model, "permittivity", opts);
  
  // alternative way to create internal submodels:
  //
  // PermittivityModel* mod = PhysicalModel::create("permittivity", opts);
  // add_submodel("permittivity", mod)

  // NOTE: all submodels are initialized automatically before calling do_init()
}
