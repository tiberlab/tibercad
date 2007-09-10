#include "HeatModel.h"
#include "Material.h"
#include "LatticeThermalConductivity.h"
#include "SimulationInterface.h"


HeatModel::HeatModel() :
  kappa(NULL),
  kappa_carrier(NULL),
  _eTEpower(0),
  _hTEpower(0),
  _elem(NULL),
  _dd_simul(NULL),
  _lattice_thermal_conductivity(0),
  _electrons_thermal_conductivity(0),
  _holes_thermal_conductivity(0)
{
}
	
 
 

HeatModel::~HeatModel()
{

  PhysicalModelInterface::destroy(kappa);
  PhysicalModelInterface::destroy(kappa_carrier);

}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::do_init()
{

  //Read Model and return a pointer to drift diffusion simulation eventually

  model_opt.joule_effect = get_options().get_option("Joule_Effect", false);

  model_opt.peltier_thomson_effect = get_options().get_option("Peltier_Thomson_Effect", false);

  model_opt.particle_thermal_conductivity = get_options().get_option("Particle_thermal_conductivity", false);

  if (model_opt.joule_effect                  |
      model_opt.peltier_thomson_effect        |  
      model_opt.particle_thermal_conductivity    )
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



  //Particle thermal conductivity

   PhysicalModelInterface::destroy(kappa_carrier);
   
   it = get_options().submodels_begin("Particle_thermal_conductivity");
   end = get_options().submodels_end("particle_thermal_conductivity");

   if (it != end)
   {

     kappa_carrier = dynamic_cast<ParticleThermalConductivity*>(
	     PhysicalModelInterface::create("particle_thermal_conductivity", it->second));
     
     if (kappa_carrier == NULL)
       throw InitFailedException("Could not create particle thermal conductivity power model");
     
     

     kappa_carrier->init();  
     
     kappa_carrier->set_material(get_material());

    }

    //Getting IDS

    //std::vector< std::map< ID, double > >  dd_solution;

    //const std::set< ID > ids;
    //ids.set(_dd_simul-> convert_variable_name_to_id("ElPotential");
    //ids.set(_dd_simul-> convert_variable_name_to_id("HoPotential");

   
    
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

  update_particle_thermal_conductivity();
  
  update_thermoelectric_powers();
    
      
}

void HeatModel::update_particle_thermal_conductivity()
{
  if (model_opt.particle_thermal_conductivity)
  {
    if (kappa_carrier != NULL)
    {
      //Insert phase    
      kappa_carrier->set_temperature(_temperature);
      
      double sigma_e = _dd_simul->get_electron_conducibility(_elem);
      
      double sigma_h = _dd_simul->get_hole_conducibility(_elem);
      
      kappa_carrier->set_electrons_conducibility(sigma_e);
      
      kappa_carrier->set_holes_conducibility(sigma_h);
      
      //Update phase
      kappa_carrier->re_init(); 
      
      //Getting result phase
      kappa_carrier->get_electrons_thermal_conductivity(_electrons_thermal_conductivity);
      
      kappa_carrier->get_holes_thermal_conductivity(_holes_thermal_conductivity);
    }
    
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

void HeatModel::update_thermoelectric_powers()
{

  if (model_opt.peltier_thomson_effect)

  {
    _eTEpower = _dd_simul->get_electrons_thermoelectric_power(_elem);
    
    _hTEpower = _dd_simul->get_holes_thermoelectric_power(_elem);

  }

}





void HeatModel::get_dd_solution( std::vector<Point> g_point,
                                std::vector<DriftDiffusion::Solution>& potentials,
                                std::vector<DriftDiffusion::Currents>& currents)
{


  if (_dd_simul != NULL & _dd_simul->is_solved())
  {

  

    // _dd_simul->get_solution (_elem,g_point,ids,dd_solution)                                       )



    

    _dd_simul->get_solution(_elem,g_point,potentials);  
    
    _dd_simul->get_solution(_elem,g_point,currents);

 

  }

  else

  {
    potentials.resize( _elem->n_nodes()  ); 
    currents.resize( _elem->n_nodes() ); 
    for (unsigned int n = 0; n < _elem->n_nodes(); n++)
      { 

	potentials[n] = (0.0);
	  currents[n] = (0.0);
	 
      }

   
  }

}
