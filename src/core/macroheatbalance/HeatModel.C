#include "HeatModel.h"
#include "Material.h"
#include "LatticeThermalConductivity.h"
#include "SimulationInterface.h"


HeatModel::HeatModel() :
  kappa(NULL),
  _eTEpower(0),
  _hTEpower(0),
  _elem(NULL),
  _dd_simul(NULL),
  _lattice_thermal_conductivity(0)
{
}
	
 
 

HeatModel::~HeatModel()
{

  PhysicalModelInterface::destroy(kappa);


}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::do_init()
{

  //  Read Model and return a pointer to drift diffusion simulation eventually

  model_opt.joule_effect = get_options().get_option("Joule_Effect", false);

  model_opt.peltier_thomson_effect = get_options().get_option("Peltier_Thomson_Effect", false);


  if (model_opt.joule_effect | model_opt.peltier_thomson_effect)
    {
      
     std::string drift_diffusion_simulation = get_options().get_option("current_simulation", "no_current");

     _dd_simul = dynamic_cast< DriftDiffusion* > ( SimulationInterface::find_simulation(drift_diffusion_simulation ) );
 
      if (_dd_simul == NULL)
      throw InitFailedException("Unknown current model" +  drift_diffusion_simulation );
    }

  
  
  //Read subModels 

  //Lattice thermal conductivuty

  PhysicalModelInterface::destroy(kappa);


  ModelOptions::const_submodel_iterator it,end;
   it = get_options().submodels_begin("Lattice_thermal_condictivity");
   end = get_options().submodels_end("Lattice_thermal_condictivity");


   if (it != end)
    {

   kappa =dynamic_cast<LatticeThermalConductivity*>(
        PhysicalModelInterface::create("lat_therm_cond_" +
		get_material()->get_structure(), it->second)); 
  
   if (kappa == NULL)
      throw InitFailedException("Could not create lattice thermal conductivity model");

 
   }
   else
    {
   kappa = dynamic_cast<LatticeThermalConductivity*>(
        PhysicalModelInterface::create("lat_therm_cond_" +
        get_material()->get_structure()));

    }
 

  kappa->set_temperature(SimulationOptions::temperature);
    
  kappa->set_material(get_material());

  kappa->init();




}




//==========================================================================//
void HeatModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const HeatModel* matA = dynamic_cast< const HeatModel*> (comp_A);

  const HeatModel* matB = dynamic_cast< const HeatModel*> (comp_B);

 
  kappa->build_alloy(matA->kappa, matB->kappa, xa);

 
}



void HeatModel::re_init()
{
  
  update_lattice_thermal_conductivity();
  

  if (model_opt.peltier_thomson_effect)
  {
    update_electrons_thermoelectric_powers();
    
    update_holes_thermoelectric_powers();
    
  }
  
 
      
}


void HeatModel::update_lattice_thermal_conductivity()
{
  //Insert phase
  kappa->set_temperature(_temperature);

  //Update phase
  kappa->update_tensor();

  //Getting result phase
  kappa->get_conductivity(_lattice_thermal_conductivity);

}

void HeatModel::update_electrons_thermoelectric_powers()
{
 
 _eTEpower = _dd_simul->get_electrons_thermoelectric_power(_elem);
  
}



void HeatModel::update_holes_thermoelectric_powers()
{

  
  _hTEpower = _dd_simul->get_holes_thermoelectric_power(_elem);
  
 
}

void HeatModel::get_dd_solution( std::vector<Point> g_point,
                                std::vector<DriftDiffusion::Solution>& potentials,
                                std::vector<DriftDiffusion::Currents>& currents)
{


  if (_dd_simul != NULL)
  {
    _dd_simul->get_solution(_elem,g_point,potentials);  
    
    _dd_simul->get_solution(_elem,g_point,currents);

 

  }
  else
  {
       for (unsigned int n = 0; n < _elem->n_nodes(); n++)
       {   
	 potentials[n] = 0.0;
	 currents[n] = 0.0;
   
       }

   
  }

}
