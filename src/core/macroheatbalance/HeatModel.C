#include "HeatModel.h"
#include "Material.h"
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

  PhysicalModelInterface::destroy(kappa);
  
  PhysicalModelInterface::destroy(thermoelectric_power);

  const ModelOptions& opt =  get_options ();

  kappa =dynamic_cast<LatticeThermalConductivity*>
    (  PhysicalModelInterface::create("lat_therm_cond_" + get_material()->get_structure(),  opt  )  );
  
  
  thermoelectric_power =dynamic_cast<ThermoelectricPower*>
    (  PhysicalModelInterface::create("thermoelectric_power" ,  opt  ))  ;



  kappa->set_material(get_material());

  kappa->init();


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

