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

void
ThermalModel::prepare_submodels(void)
{
  
  ModelOptions opts;
  opts.set_option("type","constant");
  PhysicalModelInterface* dummy;
  create_submodel(dummy, "thermal_conductivity", opts);

}

void
ThermalModel::do_init(void)
{

   PhysicalModelInterface::SubmodelIterator  it;
  //Lattice Thermal Conductivity
  it = submodels_begin("thermal_conductivity");  
  _tcm = dynamic_cast<ThermalConductivityModel*> ((*it).second);
  _kappa = _tcm->get_thermal_conductivity();
  //---------------------------------------
  ModelOptions::submodel_iterator
    it_hs(get_options().submodels_begin("heat_source"));
  ModelOptions::submodel_iterator
    end_hs(get_options().submodels_end("heat_source"));

  //Heat Source
  it = submodels_begin("heat_source");
  const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("heat_source"));
  for ( ; it != it_end ; ++it)
    _hsm.push_back(dynamic_cast<HeatSourceModel*> ((*it).second));
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
  os <<"  Kxx: "<<_kappa(0,0)<<" W/cm K\n";
  os <<"  Kyy: "<<_kappa(1,1)<<" W/cm K\n";
  os <<"  Kzz: "<<_kappa(2,2)<<" W/cm K";
  Messages::info(os.str());
  Messages::newline();
}
