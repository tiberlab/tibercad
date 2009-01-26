#include "PoissonModel.h"
#include "Material.h"
#include "ChargeDensityModel.h"
#include "Database.h"
#include "getpot.h" 



PoissonModel::PoissonModel(void) :
  _epsilon(0),
  _pyropolarization(0),
  _piezopolarization(0),
  _piezo_sim(NULL),
  _chd_sim(NULL),
  // chd_model(NULL),
  dielectric_model(NULL)
{
}
	
 

PoissonModel::~PoissonModel(void)
{

  //PhysicalModelInterface::destroy(chd_model);

PhysicalModelInterface::destroy(dielectric_model);

}

//==========================================================================//

PhysicalModelInterface* PoissonModel::create_new (void) const
{
   return new PoissonModel;
}

//==========================================================================//
void PoissonModel::do_init()
{
    //Pyropolarization
    model_opt.pyro_pol = get_options().get_option("Pyropolarization", false);

    if (model_opt.pyro_pol)
    {

       const Material* mat = get_material();
       GetPot data((mat->get_database()).get_data_file());

       _pyropolarization(1) = data("Px",0.0);  
       _pyropolarization(2) = data("Py",0.0); 
       _pyropolarization(3) = data("Pz",0.0); 
        
        const RotatedCrystal&   cr = mat->get_rotated_crystal ();
 
        _pyropolarization =cr.RotMatrix * _pyropolarization;
        
     }

   
    //Piezopolarization
    model_opt.piezo_pol = get_options().get_option("Piezopolarization", false); 
    if (model_opt.piezo_pol)
    {

      std::string piezo_sim_name = get_options().get_option("piezoelectricity_simulation", "no_sim");
      _piezo_sim =  dynamic_cast<Macrostrain* > ( SimulationInterface::find_simulation(piezo_sim_name));
      
      if (_piezo_sim == NULL) 
	throw InitFailedException("Unknown strain simulation: " +  piezo_sim_name );
      

    }


   //Density Charge Model

    model_opt.chd_sim = get_options().get_option("Charge_density_simulation", true); 

    //if (model_opt.chd_sim)
    //{
      std::string chd_sim_name = get_options().get_option("charge_density_simulation", "no_sim");
      _chd_sim =  SimulationInterface::find_simulation(chd_sim_name);
      
      if (_chd_sim == NULL) 
      	throw InitFailedException("Unknown charge_density simulation: " +  chd_sim_name );
      else
        charge_id = _chd_sim->get_variable_id("charge_density"); 
      std::cout<<charge_id<<std::endl;
      
      ModelOptions::const_submodel_iterator it,end; 

    //    PhysicalModelInterface::destroy(chd_model); */



/*    it = get_options().submodels_begin("Charge_Density_Model"); */
/*    end = get_options().submodels_end("Charge_Density_Model"); */


/*    if (it != end) */
/*     { */

/*       chd_model = dynamic_cast<ChargeDensityModel*>( */
/*         PhysicalModelInterface::create("charge_density_model", it->second));  */
  
/*       if (chd_model == NULL) */
/* 	throw InitFailedException("Could not create charge density model"); */
      
      
/*     } */
/*    else */
/*    { */
    
/*      //The Simulation name is given from Physics Section (temporany solution) */
 
/*      chd_model = dynamic_cast<ChargeDensityModel*>( */
/*      PhysicalModelInterface::create("charge_density_model"),get_options());  */
     

/*    } */

/*   chd_model->set_material(get_material()); */

/*   chd_model->init(); */


  //Dielectric Model
  PhysicalModelInterface::destroy(dielectric_model);

   it = get_options().submodels_begin("Dielectric_Model");
   end = get_options().submodels_end("Dielectric_Model");


   if (it != end)
    {

      dielectric_model = dynamic_cast< DielectricModel*>(
        PhysicalModelInterface::create("dielectric_model", it->second)); 
  
      if ( dielectric_model == NULL)
	throw InitFailedException("Could not create dielectric model");
      
      
    }
   else
   {
     
     dielectric_model = dynamic_cast<DielectricModel*>(
	 PhysicalModelInterface::create("dielectric_model")); 
     
  
   }

   dielectric_model ->set_material(get_material());
   
   dielectric_model->init();
   
   _epsilon = dielectric_model->get_dielectric_constant();



}




//==========================================================================//
void PoissonModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

   const PoissonModel* matA = dynamic_cast< const PoissonModel*> (comp_A);
   const PoissonModel* matB = dynamic_cast< const PoissonModel*> (comp_B);
  
   PhysicalModelInterface::destroy(dielectric_model);
   dielectric_model = create_submodel_alloy(matA->dielectric_model,
       matB->dielectric_model,xa);

   alloy( _pyropolarization(1),matA->_pyropolarization(1),matB->_pyropolarization(1),xa);
   alloy( _pyropolarization(2),matA->_pyropolarization(2),matB->_pyropolarization(2),xa);
   alloy( _pyropolarization(3),matA->_pyropolarization(3),matB->_pyropolarization(3),xa);
 
}



void  PoissonModel::re_init()
{
  
  // update_charge_density();

  update_built_in_polarization();

  //  update_dielectric_constant();      

}

void  PoissonModel::update_charge_density()
{ 

 
  // chd_model->set_element(_elem);

  // chd_model->re_init();

  //_charge_density = chd_model->get_charge_density();

}

//void  PoissonModel::update_dielectric_constant()
//{ 

 
//  chd_model->set_element(_elem);

//  dielectric_model->re_init();

//  charge_density = chd_model->get_charge_density();


//}

void  PoissonModel::update_built_in_polarization()
{

 if (model_opt.piezo_pol)
  { 
    SimulationEnvironment& se = _piezo_sim->get_environment();
    
    if  (se.contains_element(_elem))
    {
      
      _piezopolarization =  _piezo_sim->get_piezopolarization(_elem);

      
    }
    else
    {
      _piezopolarization=0;
    }
  }
  else
  {
    _piezopolarization=0;
  }


}

 
SimulationEnvironment& PoissonModel::get_piezo_environment()
{

  if (model_opt.piezo_pol)
  {return _piezo_sim->get_environment();}
  // else
  // {return NULL;}

  

}

  
void
PoissonModel::get_charge_density(const std::vector<Point> q_point, std::vector<double>& charge_density)
{

  
      
  if (_chd_sim->get_solution(_elem,q_point,charge_id, charge_density))
  {  
      
  }
  else
  {
     charge_density.resize(_elem->n_nodes());
     charge_density.resize(q_point.size() );
     charge_density.clear();
   }
  
      //  const Material* mat = get_material();
      //charge_density.resize(_elem->n_nodes());
   
      //for (unsigned int n = 0; n < _elem->n_nodes(); n++)
      // {
      // charge_density[n] = mat->get_net_doping_density()*Constants::e;
      // }    
  
}
