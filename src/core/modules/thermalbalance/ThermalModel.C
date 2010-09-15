// $Id$

#include "ThermalModel.h"
#include "SimulationOptions.h"
#include "Material.h"
#include "HeatSourceModel.h"
#include "ThermalConductivityModel.h"


using namespace std;



ThermalModel::ThermalModel(const ModelOptions& options)
  : PhysicalModel(options),
    _htm(NULL),
    _tcm(NULL),
    _kappa(0),
    _vg(1),
    _cg(1),
    _heat_source(0)
{
  
}

void
ThermalModel::create_submodels(void)
{

  

  //Heat Transport Default
  if (!get_options().has_submodel("heat_transport"))
  {
    ModelOptions opts;
    opts.set_option("type","fourier");
    get_options().add_submodel("heat_transport",opts);  
  }

   //Thermal Conductivity Default
   if (!get_options().has_submodel("thermal_conductivity"))
   {
     ModelOptions opts;
     opts.set_option("type","constant");
     get_options().add_submodel("thermal_conductivity",opts);
   }

}


void
ThermalModel::do_init(void) 
{

  PhysicalModelInterface::SubmodelIterator  it;
 
  //Heat Transport Model
  it = submodels_begin("heat_transport");
  _htm = dynamic_cast<HeatTransportModel*> ((*it).second);
  
 
  ModelOptions& opts = _htm->get_options();

  //Lattice Thermal Conductivity
  it = submodels_begin("thermal_conductivity");
  _tcm = dynamic_cast<ThermalConductivityModel*> ((*it).second);

 
  //Heat Source
  it = submodels_begin("heat_source");
  const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("heat_source"));
  for ( ; it != it_end ; ++it)
    _hsm.push_back(dynamic_cast<HeatSourceModel*> ((*it).second));
 

  //If we do gray we have to get the sound velocity and the relaxation time
  _kappa = _tcm->get_thermal_conductivity();
  if (_htm->get_type() == HeatTransportModel::Gray)
  {
    //----------------------------------
    Database& db = get_database();
    //Sound Velocity
    db.set_section("sound_velocity/constant");
    _vg = db.get("vg",0.0, true);
    get_parameter("vg", _vg);

    //Heat Capacity
    db.set_section("heat_capacity/constant");
    _cg = db.get("C",0.0, true);
    get_parameter("C", _cg);
   
    
    _tg = 3.0 * _kappa(0,0) / _vg / _vg / _cg;
  }

}


void
ThermalModel::do_init_alloy(const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{



   const ThermalModel* modA = dynamic_cast<const ThermalModel*>(comp_A);
   const ThermalModel* modB = dynamic_cast<const ThermalModel*>(comp_B);

   _vg = alloy(modA->_vg, modB->_vg, xa);
   _cg = alloy(modA->_cg, modB->_cg, xa);
  
}

 





void
ThermalModel::do_print_info(void)
{
  if  (SimulationOptions::verbose() > 2)
  {

    cout<<endl;
    cout<<"kx : "   <<_kappa(0,0) <<" W/cm/K"<<endl;
    cout<<"kz : "   <<_kappa(2,2) <<" W/cm/K"<<endl;
    if (_htm->get_type() == HeatTransportModel::Gray)
    {
      cout<<"tg: "  <<_tg * 1e12   <<" ps"<<endl;
      cout<<"Lg: "<<_tg * _vg * 1e7  <<" nm"<<endl;
      cout<<"vg: " <<_vg          <<" cm/s"<<endl;
      cout<<"C : "  <<_cg          <<" j/K/cm3"<<endl;
      cout<<endl;
    }
  }

  
}

void
ThermalModel::read_database(void)
{



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

  _tg = 3.0 * _kappa(0,0) / _vg / _vg / _cg;
}

