// $Id$

#include "MDModel.h"
#include "Material.h"
using namespace std;


MDModel::MDModel(const ModelOptions& options)
  : PhysicalModel(options)
{

}


MDModel*
MDModel::create(const Material* mat, const ModelOptions& options)
{


    return PhysicalModelInterface::create<MDModel>(_create,_destroy, mat, options);



  //string type("default");
  //options.get_option("type", type);

  //PhysicalModelInterface* pm = NULL;
  //pm = PhysicalModelInterface::create(_create, _destroy, options);


  //if (type == "default")
    // we create the default model from explicit creation method
  //  pm = PhysicalModelInterface::create(_create, _destroy, options);
  //else
  //{
  //  type = "bulk_" + type;
  //  pm = PhysicalModelInterface::create(type, options);
 // }

  //return dynamic_cast<MDModel*>(pm);
}



void
MDModel::do_init(void)
{
  //SubmodelIterator it = submodels_begin("charge_density");
  //if (it != submodels_end("charge_density"))
  //  _charge_density = dynamic_cast<MyChargeDensityModel*>(it->second);
  
  //it = submodels_begin("polarization");
  //const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("polarization"));
  //for ( ; it != it_end ; ++it)
  //  _pm.push_back(dynamic_cast<PolarizationModel*> ((*it).second));

  //it = submodels_begin("permittivity");
  //if (it != submodels_end("permittivity"))
   // _permittivity_model = dynamic_cast<PermittivityModel*>(it->second);

  //_permittivity = _permittivity_model->get_permittivity();
}

void
MDModel::do_calculate(void)
{
 

}

void
MDModel::create_submodels(void)
{
  
  //Thermal Conductivity Default
//  if (!get_options().has_submodel("permittivity"))
 // {
 //   ModelOptions opts;
 //   opts.set_option("type","constant");
 //   get_options().add_submodel("permittivity",opts);
 // }
  
}
