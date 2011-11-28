// $Id: ThermalModel.C 2457 2011-03-06 23:52:12Z gromano $

#include "ThermalModel.h"
#include "SimulationOptions.h"
#include "Material.h"
#include "Database.h"
#include "Messages.h"
#include "HeatSourceModel.h"
#include "ThermalConductivityModel.h"


using namespace std;



ThermalModel::ThermalModel(const ModelOptions& options)
  : PhysicalModel(options),
    _tcm(NULL),
    _kappa(0),
    _heat_source(0)
{

}


ThermalModel*
ThermalModel::create(const Material* mat, const ModelOptions& options)
{
  return PhysicalModelInterface::create<ThermalModel>(_create,_destroy, mat, options);
}

void
ThermalModel::prepare_submodels(void)
{
  
  ModelOptions opts;
  opts.set_option("type","constant");
  create_submodel(_tcm, "thermal_conductivity", opts);
  create_submodels(_hsm, "heat_source");

}

void
ThermalModel::do_init(void)
{
  _kappa = _tcm->get_thermal_conductivity();
}



void
ThermalModel::calculate(const Elem* elem, const Point& point)
{
  //Heat Source
  _heat_source = 0.0;
 
  for (ID n = 0 ; n <_hsm.size() ; n++)
  {
    _hsm[n]->calculate(elem,point);
    _heat_source +=  _hsm[n]->get_heat_source();
  }

  _tcm->calculate(elem,point);
  _kappa = _tcm->get_thermal_conductivity();

}


    //! Print some useful information
void 
ThermalModel::do_print_info(void)
{
  Messages::info("Thermal conductivity:");
  ostringstream os;
  os <<"  Kxx: "<<_kappa(0,0)<<" W/(m K)\n";
  os <<"  Kyy: "<<_kappa(1,1)<<" W/(m K)\n";
  os <<"  Kzz: "<<_kappa(2,2)<<" W/(m K)";
  Messages::info(os.str());
  Messages::newline();
}
