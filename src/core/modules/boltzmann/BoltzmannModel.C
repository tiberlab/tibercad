// $Id: BoltzmannModel.C 2457 2011-03-06 23:52:12Z gromano $

#include "BoltzmannModel.h"
#include "SimulationOptions.h"
#include "Material.h"
#include "Database.h"
#include "HeatSourceModel.h"
#include "ThermalConductivityModel.h"


using namespace std;



BoltzmannModel::BoltzmannModel(const ModelOptions& options)
  : PhysicalModel(options),
    _htm(NULL),
    _tcm(NULL),
    _kappa(0),
    _vg(1),
    _cg(1),
    _heat_source(0)
{

}

BoltzmannModel*
BoltzmannModel::create(const Material* mat, const ModelOptions& options)
{
  return PhysicalModelInterface::create<BoltzmannModel>(_create, _destroy, mat, options);
}


void
BoltzmannModel::prepare_submodels(void)
{
  ModelOptions opts;

  //Heat Transport Default
  opts.set_option("type", "fourier");
  create_submodel(_htm, "heat_transport", opts);

  
  //Thermal Conductivity Default
  opts.set_option("type","constant");
  create_submodel(_tcm, "thermal_conductivity", opts);


  create_submodels(_hsm, "heat_source");

}

void
BoltzmannModel::do_init(void)
{

  //If we do gray we have to get the sound velocity and the relaxation time
  _kappa = _tcm->get_thermal_conductivity();
  if (_htm->get_type() == HeatTransportModel::Gray)
  {
    //----------------------------------
    const Database& db = get_database();
    //Sound Velocity
    db.set_section("sound_velocity/constant");
    _vg = db.get("vg",0.0, false);
    get_parameter("vg", _vg);

    //Heat Capacity
    db.set_section("heat_capacity/constant");
    _cg = db.get("C",0.0, false);
    get_parameter("C", _cg);

    //Avarage Kappa
    double kappam = (_kappa(0,0) + _kappa(0,0)  + _kappa(0,0))/3.0;

    double mfp(0);
    get_parameter("mfp", mfp); //units: meters

    if (mfp != 0.0)
      _vg = kappam * 3.0/_cg/mfp;

    _tg = 3.0 * _kappa(0,0) / _vg / _vg / _cg;
 
  }

}


void
BoltzmannModel::do_init_alloy(const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{

   const BoltzmannModel* modA = dynamic_cast<const BoltzmannModel*>(comp_A);
   const BoltzmannModel* modB = dynamic_cast<const BoltzmannModel*>(comp_B);

   _vg = alloy(modA->_vg, modB->_vg, xa);
   _cg = alloy(modA->_cg, modB->_cg, xa);

}



void
BoltzmannModel::do_print_info(void)
{
  if  (SimulationOptions::verbose() > 2)
  {

    cout<<endl;
    cout<<"kx : "   <<_kappa(0,0) <<" W/m/K"<<endl;
    cout<<"kz : "   <<_kappa(2,2) <<" W/m/K"<<endl;
    if (_htm->get_type() == HeatTransportModel::Gray)
    {
      cout<<"tg: "  <<_tg * 1e12   <<" ps"<<endl;
      cout<<"Lg: "<<_tg * _vg<<" m"<<endl;
      cout<<"vg: " <<_vg          <<" m/s"<<endl;
      cout<<"C : "  <<_cg          <<" j/K/m3"<<endl;
      cout<<endl;
    }
  }


}


void
BoltzmannModel::calculate(const Elem* elem, const Point& point)
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

