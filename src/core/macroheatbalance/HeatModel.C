#include "HeatModel.h"
#include "Material.h"
#include "LatticeThermalConductivity.h"
HeatModel::HeatModel() 
{
	
  kappa = NULL;
  thermoelectric_power = NULL;
 
}


HeatModel::~HeatModel()
{

  PhysicalModelInterface::destroy(kappa);
  
  PhysicalModelInterface::destroy(thermoelectric_power);

}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::do_init()
{

  const ModelOptions& opt =  get_options ();
  //Kappa

  PhysicalModelInterface::destroy(kappa);

  ModelOptions::const_submodel_iterator  it = get_options().submodels_begin("Lattice_thermal_conductivity");

  const ModelOptions& kappa_option = (it->second);

  kappa =dynamic_cast<LatticeThermalConductivity*>
    (  PhysicalModelInterface::create("lat_therm_cond_" + get_material()->get_structure(),  kappa_option  )  );

  kappa->temperature = SimulationOptions::temperature;
    
  kappa->set_material(get_material());

  kappa->init();


  //Thermoelectric_power

  PhysicalModelInterface::destroy(thermoelectric_power);


  thermoelectric_power =dynamic_cast<ThermoelectricPower*>
    (  PhysicalModelInterface::create("thermoelectric_power" ,  opt  ))  ;

  thermoelectric_power->set_material(get_material());
  
  thermoelectric_power->init();  //The method init is non virtual in PMI which calls read_database and do_init that are virtual


}




//==========================================================================//
void HeatModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const HeatModel* matA = dynamic_cast< const HeatModel*> (comp_A);

  const HeatModel* matB = dynamic_cast< const HeatModel*> (comp_B);

 

  kappa->build_alloy(matA->kappa, matB->kappa, xa);

//   thermoelectric_power_n->build_alloy(matA->kappa, matB->kappa, xa);
//   thermoelectric_power_p->build_alloy(matA->kappa, matB->kappa, xa);
  

 
}

