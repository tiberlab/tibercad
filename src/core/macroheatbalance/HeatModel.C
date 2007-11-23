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

    
     _dd_simul = SimulationInterface::find_simulation(drift_diffusion_simulation);
     
    
     
     if (_dd_simul == NULL)
       throw InitFailedException("Unknown drift diffusion simulation" +  drift_diffusion_simulation );

    }
  

   model_opt.excitons = get_options().get_option("Excitons_dissipation", false);

   if (model_opt.excitons)
    {

     std::string excitons_name_simulation = get_options().get_option("excitons_simulation", "no_excitons");

     _ex_simul = SimulationInterface::find_simulation(excitons_name_simulation);

      if ( _ex_simul == NULL)
          throw InitFailedException("Unknown excitons simulation" + excitons_name_simulation);

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

     //std::cout<<
     //get_material()->get_name()<<std::endl;
  
     
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

    if (model_opt.joule_effect) 
    {
     dd_ID_je.insert(_dd_simul->get_variable_id("QFermi_e"));
     dd_ID_je.insert(_dd_simul->get_variable_id("QFermi_h"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jn_x"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jn_y"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jn_z"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jp_x"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jp_y"));
     dd_ID_je.insert(_dd_simul->get_variable_id("Jp_z"));
  
     ID_je.resize(8);
     ID_je[QFERMIE]=_dd_simul->get_variable_id("QFermi_e"); 
     ID_je[QFERMIH]=_dd_simul->get_variable_id("QFermi_h");
     ID_je[JNX]=_dd_simul->get_variable_id("Jn_x");
     ID_je[JNY]=_dd_simul->get_variable_id("Jn_y");
     ID_je[JNZ]=_dd_simul->get_variable_id("Jn_z");
     ID_je[JPX]=_dd_simul->get_variable_id("Jp_x");
     ID_je[JPY]=_dd_simul->get_variable_id("Jp_y");
     ID_je[JPZ]=_dd_simul->get_variable_id("Jp_z");
   } 
 
   if (model_opt.particle_thermal_conductivity) 
   {   

   dd_ID_kpart.insert(_dd_simul->get_variable_id("CondE")); 
   dd_ID_kpart.insert(_dd_simul->get_variable_id("CondH")); 

   ID_kpart.resize(2);
   ID_kpart[CONDE] = (_dd_simul->get_variable_id("CondE")); 
   ID_kpart[CONDH] = (_dd_simul->get_variable_id("CondH")); 
   } 

   if  (model_opt.peltier_thomson_effect) 
   { 

 
   dd_ID_TEpower.insert(_dd_simul->get_variable_id("TEpowerE")); 
   dd_ID_TEpower.insert(_dd_simul->get_variable_id("TEpowerH")); 

   ID_TEpower.resize(2);
   ID_TEpower[TEPOWERE] = (_dd_simul->get_variable_id("TEpowerE")); 
   ID_TEpower[TEPOWERH] = (_dd_simul->get_variable_id("TEpowerH")); 
      
   }  

   if (model_opt.excitons) 
   {
    
     ex_set_ID.insert(_ex_simul->get_variable_id("J_x"));
     ex_set_ID.insert(_ex_simul->get_variable_id("J_y"));
     ex_set_ID.insert(_ex_simul->get_variable_id("J_z"));
     ex_set_ID.insert(_ex_simul->get_variable_id("chemPot"));
     ex_set_ID.insert(_ex_simul->get_variable_id("Rad_power"));
  
     ID_ex.resize(5);
     ID_ex[JEX_X]=_ex_simul->get_variable_id("J_x"); 
     ID_ex[JEX_Y]=_ex_simul->get_variable_id("J_y");
     ID_ex[JEX_Z]=_ex_simul->get_variable_id("J_z");
     ID_ex[EX_POTENTIAL]=_ex_simul->get_variable_id("chemPot");
     ID_ex[RADPOWER]=_ex_simul->get_variable_id("Rad_power");
   }

 
   
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

      std::vector< std::map< ID, double > >  dd_sol_kpart;
      std::vector<Point> centroid;
      centroid[0]=_elem->centroid();
    
      _dd_simul->get_solution(_elem,centroid,dd_ID_kpart,dd_sol_kpart); 

      double sigma_e =  dd_sol_kpart[0].find(ID_kpart[CONDE])->second;
  
      double sigma_h =  dd_sol_kpart[0].find(ID_kpart[CONDH])->second;
      
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

  std::vector< std::map< ID, double > >  dd_sol_te;
  std::vector<Point> centroid;
  centroid[0]=_elem->centroid();
  
  _dd_simul->get_solution(_elem,centroid,dd_ID_TEpower,dd_sol_te); 

  _eTEpower =  dd_sol_te[0].find(ID_TEpower[TEPOWERE])->second; 
      
  _hTEpower =  dd_sol_te[0].find(ID_TEpower[TEPOWERH])->second;   


  }

}

 
void HeatModel::get_dd_solution_secure( std::vector<Point> g_point,
                                        std::vector<double>& QfermiE,
                                        std::vector<double>& QfermiH,
                                        std::vector<Point>& JE,
                                        std::vector<Point>& JH)
{

     std::vector< std::map< ID, double > >  dd_sol_je;

     _dd_simul->get_solution(_elem,g_point,dd_ID_je,dd_sol_je); 
         
     for (unsigned int n=0; n<g_point.size(); n++) // loop over test function
     {

       QfermiE.resize(g_point.size());
       QfermiH.resize(g_point.size());
       JE.resize(g_point.size());
       JH.resize(g_point.size());

       QfermiE[n]  = dd_sol_je[n].find(ID_je[QFERMIE])->second;
       QfermiH[n]  = dd_sol_je[n].find(ID_je[QFERMIH])->second;

       JE[n](0) = dd_sol_je[n].find(ID_je[JNX])->second;

       JE[n](1) = dd_sol_je[n].find(ID_je[JNY])->second;
       JE[n](2) = dd_sol_je[n].find(ID_je[JNZ])->second;
       JH[n](0) = dd_sol_je[n].find(ID_je[JPX])->second;
       JH[n](1) = dd_sol_je[n].find(ID_je[JPY])->second;
       JH[n](2) = dd_sol_je[n].find(ID_je[JPZ])->second;   
     }
 
  

}

void  HeatModel::get_ex_solution_secure( std::vector<Point> g_point,
			       std::vector<double>& ex_potential,
					 std::vector<Point>& J_ex,
					 std::vector<double>& radiative_power
                                     )
{

     std::vector< std::map< ID, double > >  ex_sol;

     _ex_simul->get_solution(_elem,g_point,ex_set_ID,ex_sol); 
         
     for (unsigned int n=0; n<g_point.size(); n++) // loop over test function
     {

       ex_potential.resize(g_point.size());

      
       J_ex.resize(g_point.size());

       ex_potential[n]  = Constants::e * ex_sol[n].find(ID_ex[EX_POTENTIAL])->second;
     

       J_ex[n](0) = ex_sol[n].find(ID_ex[JEX_X])->second;
       J_ex[n](1) = ex_sol[n].find(ID_ex[JEX_Y])->second;
       J_ex[n](2) = ex_sol[n].find(ID_ex[JEX_Z])->second;

       radiative_power.resize(g_point.size());
       radiative_power[n] =  ex_sol[n].find(ID_ex[RADPOWER])->second;
 

      
     }
 
  

}
